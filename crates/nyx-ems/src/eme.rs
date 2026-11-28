// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! 3D Finite-Element Eigenmode Expansion (`Study::Eme`, `"eme"`, ADR 140) along wafer `+y` (`exp(-j beta y)`).
//!
//! Structured in the three stages of ADR 9:
//! 1. [`assemble_eme`]: partitions `y in [y0, y1]` according to [`EmeGridSpec`] (`Uniform`, `Explicit`, `Composite`,
//!    `PeriodicSubgrid`), extracts 2D `(x, z)` cross-section meshes (preserving vertical-prism base triangulations on
//!    swept/extruded slabs and slicing arbitrary sloped/unstructured 3D tetrahedral meshes cleanly), evaluates
//!    domain-filtered local temperatures `T(x, y_c, z)` and potentials `V(x, y_c, z)` via [`femcore::transfer`],
//!    and deduplicates identical slices by canonical BLAKE3 digest so repeated cross-sections solve **once**.
//! 2. [`solve_eme_cascade`]: solves `M_c` guided/evanescent eigenmodes per unique slice (preserving PMLs, midplane
//!    symmetry, gyrotropic tensors, zero-thickness `sheet_2d` layers, thermo-optic `dn/dT`, and exact Ward–Pendry
//!    cylindrical bend metrics `h_y(x) = 1 + x / R_bend`), computes biorthogonal / cross-power interface overlap
//!    matrices, builds unconditionally stable interface scattering matrices `S_if` and passive propagation matrices
//!    `P_c = diag(exp(-j beta_{c,m} L_c))` (`Im(beta) <= 0`), and cascades the device via Redheffer star products
//!    (with `O(log N_rep)` binary exponentiation on `PeriodicSubgrid` nodes).
//! 3. [`eme_fields`]: emits [`DrivenPortData`] (`S`, `Z`, `Y` matrices compatible with `nyx_analysis::network` and
//!    Touchstone export), caches [`EmeData`] in `SolveResult::scalars` for microsecond zero-re-solve post-processing
//!    ([`eme_length_sweep`], [`eme_periodicity_sweep`], [`eme_mode_sweep`]), and reconstructs the 3D electric field
//!    via bounded forward/backward prefix-suffix scattering amplitudes.

use crate::ibc::{add_ibc_2d_waveguide, add_sheet2d_2d_waveguide, pec_edges_and_nodes_2d, MU0};
use crate::mixed_ops::{apply_thermo_optic_to_c3, c3_mul_vec, is_c3_slice_non_symmetric, EmPmlTable, C3};
use crate::modes::{assemble_modes_with_bcs_and_symmetry, solve_mode_pairs, EpsCell, C0};
use crate::modes_complex::{add_waveguide_complex_point_full, solve_mode_pairs_complex, ModeAssemblyC, ModesC};
use crate::problem::{pec_edges_and_nodes_2d_with_symmetry, Problem, Sheet2dBc, Study, Symmetry, WavePortBc};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::{free_dofs, geometry, quad_rule, shape};
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{c, CCoo, C};
use femcore::linalg::prom::{s_matrix_passivity_and_reciprocity, s_to_z_and_y_checked, solve_dense_complex_multi};
use femcore::mesh::SimplexMesh;
use femcore::optical::fine_rule;
use femcore::props::Pml;
use femcore::section::{local_tensor, PROP_SIGN};
use femcore::transfer::{canonical_name, SourceField, SourceMesh, TransferIndex, TransferOptions};
use nyx_physics::decode::Params;
use nyx_physics::{
    cell_fields_of, discretise, Discretisation, DrivenPortData, Error, Field, Mesh, Scalar, Solution, SolveResult,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::f64::consts::PI;
use std::time::Instant;

/// Scalar key carrying serialized [`EmeData`] on `ems` `eme` results for zero-re-solve sweeps (`ADR 140`).
pub const EME_DATA_KEY: &str = "eme cascade data";

/// Allowed flat choices for `grid_kind` on `StudyDesc` `"eme"`.
pub const EME_GRID_KIND_CHOICES: [&str; 5] = ["uniform", "explicit", "periodic_subgrid", "composite", "json"];

/// Hierarchical partition of the propagation interval `y in [y0, y1]` into translationally invariant EME segments.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum EmeGridSpec {
    /// `n_cells` equal-length segments across the assigned `y` interval.
    Uniform {
        n_cells: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        n_modes: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        bend_radius: Option<f64>,
    },
    /// Explicit sequence of segment lengths `lengths` (metres) along `+y`.
    Explicit {
        lengths: Vec<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        n_modes: Option<Vec<usize>>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        bend_radius: Option<Vec<Option<f64>>>,
        /// Optional fixed mesh-sampling weights along `+y` when physical propagation `lengths` are swept/optimized independently of the 3D mesh.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        mesh_lengths: Option<Vec<f64>>,
    },
    /// Concatenation of child subgrids along `+y`.
    Composite { subgrids: Vec<EmeGridSpec> },
    /// Periodic repetition (`num_reps >= 1`) of `unit_cell`, evaluated in `O(log num_reps)` Redheffer star products.
    PeriodicSubgrid {
        unit_cell: Box<EmeGridSpec>,
        num_reps: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        unit_length: Option<f64>,
    },
}

impl EmeGridSpec {
    /// True when every numeric length and bend radius in the specification is finite.
    pub fn is_finite(&self) -> bool {
        match self {
            Self::Uniform { bend_radius, .. } => bend_radius.is_none_or(f64::is_finite),
            Self::Explicit { lengths, bend_radius, mesh_lengths, .. } => {
                lengths.iter().all(|l| l.is_finite())
                    && bend_radius.as_ref().is_none_or(|bs| bs.iter().all(|b| b.is_none_or(f64::is_finite)))
                    && mesh_lengths.as_ref().is_none_or(|ml| ml.iter().all(|l| l.is_finite()))
            }
            Self::Composite { subgrids } => subgrids.iter().all(Self::is_finite),
            Self::PeriodicSubgrid { unit_cell, unit_length, .. } => {
                unit_cell.is_finite() && unit_length.is_none_or(f64::is_finite)
            }
        }
    }

    fn explicit_optical_length(&self) -> Option<f64> {
        match self {
            Self::Uniform { .. } => None,
            Self::Explicit { lengths, .. } => Some(lengths.iter().sum()),
            Self::Composite { subgrids } => {
                let mut total = 0.0;
                for sg in subgrids {
                    total += sg.explicit_optical_length()?;
                }
                Some(total)
            }
            Self::PeriodicSubgrid { unit_cell, num_reps, unit_length } => {
                let lu = unit_length.or_else(|| unit_cell.explicit_optical_length())?;
                Some(lu * (*num_reps as f64))
            }
        }
    }

    fn explicit_unit_mesh_length(&self) -> Option<f64> {
        match self {
            Self::Uniform { .. } => None,
            Self::Explicit { lengths, mesh_lengths, .. } => Some(mesh_lengths.as_ref().unwrap_or(lengths).iter().sum()),
            Self::Composite { subgrids } => {
                let mut total = 0.0;
                for sg in subgrids {
                    total += sg.explicit_unit_mesh_length()?;
                }
                Some(total)
            }
            Self::PeriodicSubgrid { unit_cell, unit_length, .. } => {
                unit_length.or_else(|| unit_cell.explicit_unit_mesh_length())
            }
        }
    }

    fn uniform_cell_weight(&self) -> usize {
        match self {
            Self::Uniform { n_cells, .. } => (*n_cells).max(1),
            Self::Explicit { .. } => 0,
            Self::Composite { subgrids } => subgrids.iter().map(Self::uniform_cell_weight).sum(),
            Self::PeriodicSubgrid { unit_cell, num_reps, unit_length } => {
                if unit_length.is_some() || unit_cell.explicit_unit_mesh_length().is_some() {
                    0
                } else {
                    unit_cell.uniform_cell_weight().saturating_mul((*num_reps).max(1))
                }
            }
        }
    }
}

/// Decodes an [`EmeGridSpec`] from the flat/JSON parameters of `StudyDesc` `"eme"`.
pub(crate) fn decode_eme_grid(sp: &Params) -> Result<EmeGridSpec, Error> {
    let grid_json = sp.text("grid_json")?;
    let grid_kind = sp.text("grid_kind")?;
    if !grid_json.trim().is_empty() {
        let parsed: EmeGridSpec = serde_json::from_str(&grid_json)
            .map_err(|e| Error::invalid(format!("eme study: invalid 'grid_json': {e}")))?;
        validate_grid_spec(&parsed).map_err(Error::invalid)?;
        return Ok(parsed);
    }
    let bend_raw = sp.number("bend_radius")?;
    if !bend_raw.is_finite() {
        return Err(Error::invalid("eme study: 'bend_radius' must be finite"));
    }
    let bend = (bend_raw.abs() > 0.0).then_some(bend_raw);
    let n_cells = sp.integer("n_cells")?.max(1) as usize;
    let lengths = sp.numbers("lengths")?;
    let num_reps = sp.integer("num_reps")?.max(1) as usize;

    let spec = match grid_kind.as_str() {
        "uniform" | "" => EmeGridSpec::Uniform { n_cells, n_modes: None, bend_radius: bend },
        "explicit" => {
            if lengths.is_empty() {
                return Err(Error::invalid("eme study: 'explicit' grid_kind requires non-empty 'lengths'"));
            }
            let br = bend.map(|r| vec![Some(r); lengths.len()]);
            EmeGridSpec::Explicit { lengths, n_modes: None, bend_radius: br, mesh_lengths: None }
        }
        "periodic_subgrid" => {
            let unit = if !lengths.is_empty() {
                let br = bend.map(|r| vec![Some(r); lengths.len()]);
                EmeGridSpec::Explicit { lengths, n_modes: None, bend_radius: br, mesh_lengths: None }
            } else {
                EmeGridSpec::Uniform { n_cells, n_modes: None, bend_radius: bend }
            };
            EmeGridSpec::PeriodicSubgrid { unit_cell: Box::new(unit), num_reps, unit_length: None }
        }
        "composite" | "json" => {
            return Err(Error::invalid(format!(
                "eme study: grid_kind '{grid_kind}' requires a non-empty 'grid_json' specification"
            )));
        }
        other => {
            return Err(Error::invalid(format!(
                "eme study: unknown grid_kind '{other}' (expected one of {EME_GRID_KIND_CHOICES:?})"
            )));
        }
    };
    validate_grid_spec(&spec).map_err(Error::invalid)?;
    Ok(spec)
}

fn validate_grid_spec(spec: &EmeGridSpec) -> Result<(), String> {
    match spec {
        EmeGridSpec::Uniform { n_cells, n_modes, bend_radius } => {
            if *n_cells == 0 {
                return Err("EME Uniform subgrid requires n_cells >= 1".into());
            }
            if let Some(m) = n_modes {
                if *m == 0 {
                    return Err("EME Uniform subgrid requires n_modes >= 1 when specified".into());
                }
            }
            if let Some(r) = bend_radius {
                if !r.is_finite() || *r == 0.0 {
                    return Err(format!("EME bend_radius must be finite and non-zero (got {r})"));
                }
            }
            Ok(())
        }
        EmeGridSpec::Explicit { lengths, n_modes, bend_radius, mesh_lengths } => {
            if lengths.is_empty() {
                return Err("EME Explicit subgrid requires at least one segment length".into());
            }
            for (i, &l) in lengths.iter().enumerate() {
                if !(l > 0.0 && l.is_finite()) {
                    return Err(format!("EME Explicit segment {i} length must be positive and finite (got {l})"));
                }
            }
            if let Some(ml) = mesh_lengths {
                if ml.len() != lengths.len() {
                    return Err(format!(
                        "EME Explicit mesh_lengths length ({}) does not match lengths ({})",
                        ml.len(),
                        lengths.len()
                    ));
                }
                for (i, &m) in ml.iter().enumerate() {
                    if !(m > 0.0 && m.is_finite()) {
                        return Err(format!(
                            "EME Explicit segment {i} mesh_lengths entry must be positive and finite (got {m})"
                        ));
                    }
                }
            }
            if let Some(ms) = n_modes {
                if ms.len() != lengths.len() {
                    return Err(format!(
                        "EME Explicit n_modes length ({}) does not match lengths ({})",
                        ms.len(),
                        lengths.len()
                    ));
                }
                if ms.contains(&0) {
                    return Err("EME Explicit n_modes entries must be >= 1".into());
                }
            }
            if let Some(bs) = bend_radius {
                if bs.len() != lengths.len() {
                    return Err(format!(
                        "EME Explicit bend_radius length ({}) does not match lengths ({})",
                        bs.len(),
                        lengths.len()
                    ));
                }
                for (i, b) in bs.iter().enumerate() {
                    if let Some(r) = b {
                        if !r.is_finite() || *r == 0.0 {
                            return Err(format!(
                                "EME Explicit segment {i} bend_radius must be non-zero finite (got {r})"
                            ));
                        }
                    }
                }
            }
            Ok(())
        }
        EmeGridSpec::Composite { subgrids } => {
            if subgrids.is_empty() {
                return Err("EME Composite grid requires at least one child subgrid".into());
            }
            for sg in subgrids {
                validate_grid_spec(sg)?;
            }
            Ok(())
        }
        EmeGridSpec::PeriodicSubgrid { unit_cell, num_reps, unit_length } => {
            if *num_reps == 0 {
                return Err("EME PeriodicSubgrid requires num_reps >= 1".into());
            }
            if let Some(lu) = unit_length {
                if !(*lu > 0.0 && lu.is_finite()) {
                    return Err(format!("EME PeriodicSubgrid unit_length must be positive finite (got {lu})"));
                }
            }
            validate_grid_spec(unit_cell)
        }
    }
}

/// Exact conformal radial coordinate `u(x) = R_bend * ln(1 + x / R_bend)` for a waveguide bent in the wafer `(x, y)` plane.
pub fn conformal_bend_coord(x: f64, r_bend: f64) -> Result<f64, String> {
    if !x.is_finite() || !r_bend.is_finite() || r_bend.abs() == 0.0 {
        return Err(format!(
            "conformal_bend_coord requires finite x and non-zero finite r_bend (got x={x}, r_bend={r_bend})"
        ));
    }
    let h = 1.0 + x / r_bend;
    if !h.is_finite() || h <= 1e-3 {
        return Err(format!(
            "bent waveguide coordinate x = {x:.3e} m reaches or crosses the bend center r_bend = {r_bend:.3e} m (1 + x / r_bend = {h:.3e} <= 1e-3)"
        ));
    }
    Ok(r_bend * h.ln())
}

/// Applies the exact Ward–Pendry cylindrical bend metric `J_b(x) = diag(1, h_y(x), 1)` (`h_y = 1 + x / R_bend`)
/// to a wafer-frame permittivity or permeability tensor `t`: `T_bend = J_b^-1 T J_b^-T det(J_b)`.
#[inline]
fn apply_bend_metric_wafer(t: &C3, h_y: f64) -> C3 {
    let h = c(h_y, 0.0);
    let inv_h = c(1.0 / h_y, 0.0);
    [[t[0][0] * h, t[0][1], t[0][2] * h], [t[1][0], t[1][1] * inv_h, t[1][2]], [t[2][0] * h, t[2][1], t[2][2] * h]]
}

/// Applies the inverse cylindrical bend metric `T_bend^-1 = (1 / h_y) J_b^T T^-1 J_b` to a wafer-frame inverse permeability `mui`.
#[inline]
fn apply_bend_metric_inv_wafer(mui: &C3, h_y: f64) -> C3 {
    let h = c(h_y, 0.0);
    let inv_h = c(1.0 / h_y, 0.0);
    [
        [mui[0][0] * inv_h, mui[0][1], mui[0][2] * inv_h],
        [mui[1][0], mui[1][1] * h, mui[1][2]],
        [mui[2][0] * inv_h, mui[2][1], mui[2][2] * inv_h],
    ]
}

/// Cached summary of one deduplicated 2D cross-section slice inside [`EmeData`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EmeUniqueSliceSummary {
    /// Hex-encoded BLAKE3 digest of the canonical slice geometry, boundary conditions, bend radius, and local material/temperature state.
    pub blake3_hash: String,
    pub sample_y_m: f64,
    pub n_modes: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bend_radius_m: Option<f64>,
    /// Forward (`+y`) effective indices `[re, im]`.
    pub neff_fwd: Vec<[f64; 2]>,
    /// Backward (`-y`) effective indices `[re, im]` (negated convention so passive decay is `Im <= 0`).
    pub neff_bwd: Vec<[f64; 2]>,
    /// Forward (`+y`) propagation constants `beta_fwd` (`rad/m`, passive `Im <= 0` guaranteed).
    pub beta_fwd: Vec<[f64; 2]>,
    /// Backward (`-y`) effective propagation constants `beta_bwd` (`rad/m`, passive `Im <= 0` guaranteed so `exp(-j beta_bwd L_c)` decays).
    pub beta_bwd: Vec<[f64; 2]>,
    /// Intra-slice biorthogonal self-pairings `D_{mm} = int (e_tilde_m^* x h_m) . y_hat dS`.
    pub d_self: Vec<[f64; 2]>,
    /// Dual intra-slice pairings `D_tilde_{mm} = int (e_m x h_tilde_m^*) . y_hat dS`.
    pub d_self_tilde: Vec<[f64; 2]>,
}

/// Cached cross-power and biorthogonal interface overlap matrices between ordered slice pair `(left_id, right_id)`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EmeInterfaceOverlap {
    pub left_id: usize,
    pub right_id: usize,
    pub n_left: usize,
    pub n_right: usize,
    pub d_left: Vec<[f64; 2]>,
    pub d_right_tilde: Vec<[f64; 2]>,
    /// `n_left x n_right` matrix `V_{LR}(m, n) = int (e_tilde_m^(L)* x h_n^(R)) . y_hat dS`.
    pub v_lr: Vec<Vec<[f64; 2]>>,
    /// `n_left x n_right` matrix `W_{LR}(m, n) = int (e_m^(L) x h_tilde_n^(R)*) . y_hat dS`.
    pub w_lr: Vec<Vec<[f64; 2]>>,
    /// Symmetric cross-power overlap `O_{c,c+1}(m, n) = 0.5 int (E_m^(c) x H_n^(c+1)* + E_n^(c+1)* x H_m^(c)) . y_hat dS`.
    pub cross_power_overlap: Vec<Vec<[f64; 2]>>,
    /// Full `(n_left + n_right) x (n_left + n_right)` interface scattering matrix `S_if`.
    pub s_interface: Vec<Vec<[f64; 2]>>,
}

/// Recursive cascade tree stored in [`EmeData`] preserving both flat segment indices and periodic subgrids.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum EmeCascadeNode {
    Segment {
        segment_index: usize,
        unique_slice: usize,
        length_m: f64,
        n_modes: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        bend_radius_m: Option<f64>,
    },
    Composite {
        children: Vec<EmeCascadeNode>,
    },
    Periodic {
        subgrid_index: usize,
        unit_cell: Box<EmeCascadeNode>,
        num_reps: usize,
    },
}

/// Complete serializable EME cascade cache stored on `SolveResult::scalars` under [`EME_DATA_KEY`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EmeData {
    pub freq_hz: f64,
    pub y_range_m: [f64; 2],
    pub total_segments: usize,
    pub total_periodic_subgrids: usize,
    pub unique_slices: Vec<EmeUniqueSliceSummary>,
    pub tree: EmeCascadeNode,
    pub interfaces: Vec<EmeInterfaceOverlap>,
    /// Optional explicit input waveport-to-slice-0 overlap (`left_id = usize::MAX - 1`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub input_port_overlap: Option<EmeInterfaceOverlap>,
    /// Optional explicit slice-(C-1)-to-output waveport overlap (`right_id = usize::MAX`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub output_port_overlap: Option<EmeInterfaceOverlap>,
    /// Mapping from output port index `0..P` to `(is_output_side, channel_index_on_side)`.
    pub port_channels: Vec<(bool, usize)>,
    pub ports: Vec<String>,
    pub port_z0_ohm: Vec<f64>,
}

impl EmeData {
    /// Reads [`EmeData`] from `result.scalars` when present.
    pub fn of(result: &SolveResult) -> Option<Self> {
        result.scalars.iter().find(|(k, _)| k == EME_DATA_KEY).and_then(|(_, v)| match v {
            Scalar::Text(t) => serde_json::from_str(t).ok(),
            Scalar::Num(_) => None,
        })
    }

    /// Encodes [`EmeData`] as a `(String, Scalar)` pair for `SolveResult::scalars`.
    pub fn scalar(&self) -> Result<(String, Scalar), String> {
        let txt = serde_json::to_string(self).map_err(|e| format!("failed to serialize EmeData: {e}"))?;
        Ok((EME_DATA_KEY.to_string(), Scalar::Text(txt)))
    }
}

/// One evaluated point in an analytical post-solve EME sweep (`length`, `periodicity`, or `modes`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EmeSweepPoint {
    /// Swept parameter value (`length` in metres, `num_reps`, or `n_modes`).
    pub parameter: f64,
    /// `P x P` complex terminal scattering matrix `[re, im]`.
    pub s_matrix: Vec<Vec<[f64; 2]>>,
    /// Power transmission `|S_{p, 0}|^2` from port 0 to every port `p`.
    pub power_from_port0: Vec<f64>,
    /// Maximum singular value `sigma_max(S)` (`<= 1` for passive cascades).
    pub passivity_max_sigma: f64,
    /// Reciprocity asymmetry `||S - S^T||_F / max(1, ||S||_F)`.
    pub reciprocity_rel_err: f64,
}

/// Output of [`eme_length_sweep`], [`eme_periodicity_sweep`], and [`eme_mode_sweep`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EmeSweepReport {
    pub mode: String,
    pub target_index: usize,
    pub freq_hz: f64,
    pub ports: Vec<String>,
    pub port_z0_ohm: Vec<f64>,
    pub points: Vec<EmeSweepPoint>,
    pub elapsed_s: f64,
}

// ------------------------------------------------------------------------------------------------------------- Stage 1: Slicing & Deduplication

/// Internal mesh-sampledsegment descriptor before deduplication.
#[derive(Clone, Debug)]
pub struct EmeSegmentPlan {
    pub segment_index: usize,
    pub y_mesh_lo: f64,
    pub y_mesh_hi: f64,
    pub y_sample: f64,
    pub length_m: f64,
    pub n_modes: usize,
    pub bend_radius_m: Option<f64>,
    pub unique_slice: usize,
}

/// Deduplicated 2D slice ready for Stage 2 eigenmode analysis.
#[derive(Clone, Debug)]
pub struct EmeUniqueSlice {
    pub unique_id: usize,
    pub blake3_hash: String,
    pub sample_y_m: f64,
    pub mesh_2d: SimplexMesh,
    pub max_modes: usize,
    pub bend_radius_m: Option<f64>,
    /// Quadrature-point temperatures `cell_quad_temps[ci][q]` (`quad_rule(2, 6)` or `quad_rule(2, 4)`), empty when uniform.
    pub cell_quad4_temps: Vec<Vec<f64>>,
    pub cell_quad6_temps: Vec<Vec<f64>>,
    /// Per-cell average temperature (K) for fast real-solver path when sub-cell variation is negligible.
    pub cell_avg_temps: Vec<f64>,
    pub has_subcell_temp_grad: bool,
    /// Quadrature-point DC electric fields `[Ex, Ey, Ez]` (V/m, wafer axes) for electro-optic perturbation, empty when inactive.
    pub cell_quad4_e_dc: Vec<Vec<[f64; 3]>>,
    pub cell_quad6_e_dc: Vec<Vec<[f64; 3]>>,
    pub cell_avg_e_dc: Vec<[f64; 3]>,
    /// Quadrature-point mechanical strains `[Sxx, Syy, Szz, Sxy, Syz, Sxz]` (`m/m`, solver Voigt) for photoelastic perturbation, empty when inactive.
    pub cell_quad4_strain: Vec<Vec<[f64; 6]>>,
    pub cell_quad6_strain: Vec<Vec<[f64; 6]>>,
    pub cell_avg_strain: Vec<[f64; 6]>,
    pub has_subcell_eo_pe_grad: bool,
    /// Optional P2 nodal electrostatic potential `V(x, y_c, z)` on `mesh_2d` for `sheet_2d` chemical-potential gating.
    pub slice_potential: Option<Vec<f64>>,
}

/// Extracted 2D waveport terminal patch on `y = y0` or `y = y1`.
#[derive(Clone, Debug)]
pub struct EmeWavePortPatch {
    pub bc: WavePortBc,
    pub is_output: bool,
    pub mesh_2d: SimplexMesh,
}

/// Stage 1 output of [`assemble_eme`] (`ADR 9`, `ADR 140`).
#[derive(Clone, Debug)]
pub struct EmeAssembly {
    pub mesh_3d: SimplexMesh,
    pub problem: Problem,
    pub freq_hz: f64,
    pub n_target: f64,
    pub evanescent: bool,
    pub y_range_m: [f64; 2],
    pub port_modes: usize,
    pub port_z0_ohm: f64,
    pub reconstruct_fields: bool,
    pub tree: EmeCascadeNode,
    pub total_periodic_subgrids: usize,
    pub mesh_segments: Vec<EmeSegmentPlan>,
    pub unique_slices: Vec<EmeUniqueSlice>,
    pub waveport_patches: Vec<EmeWavePortPatch>,
}

/// Stage 2 solved eigenmode set on a unique 2D slice (or waveport patch).
#[derive(Clone, Debug)]
pub struct SolvedSliceModes {
    pub mesh_2d: SimplexMesh,
    pub order: usize,
    pub freq_hz: f64,
    pub bend_radius_m: Option<f64>,
    pub pml_cells: Vec<bool>,
    pub mui_phys_wafer: Vec<C3>,
    pub neff_fwd: Vec<C>,
    pub neff_bwd: Vec<C>,
    pub beta_fwd: Vec<C>,
    pub beta_bwd: Vec<C>,
    pub et: Vec<Vec<C>>,
    pub phi: Vec<Vec<C>>,
    pub et_left: Vec<Vec<C>>,
    pub phi_left: Vec<Vec<C>>,
    pub d_self: Vec<C>,
    pub d_self_tilde: Vec<C>,
    pub residuals: Vec<f64>,
    pub restarts: usize,
    pub converged: bool,
    pub is_non_symmetric: bool,
}

/// Stage 2 output of [`solve_eme_cascade`] (`ADR 9`, `ADR 140`).
#[derive(Clone, Debug)]
pub struct EmeCascade {
    pub slice_modes: Vec<SolvedSliceModes>,
    pub input_port_modes: Option<SolvedSliceModes>,
    pub output_port_modes: Option<SolvedSliceModes>,
    pub data: EmeData,
    /// Full terminal `P x P` scattering matrix `S` at the nominal design point.
    pub s_terminal: Vec<Vec<C>>,
    pub max_residual: f64,
    pub total_restarts: usize,
    pub all_converged: bool,
}

/// Stage 1: partitions `y in [y0, y1]`, extracts 2D cross-section meshes and domain-filtered thermal/electrostatic
/// profiles, and deduplicates identical slices with BLAKE3.
pub fn assemble_eme(mesh: &Mesh, problem: &Problem) -> Result<EmeAssembly, String> {
    if mesh.dim() != 3 {
        return Err(format!("the study 'eme' needs a 3D mesh, this one is {}D", mesh.dim()));
    }
    let Study::Eme {
        freq,
        n_target,
        evanescent,
        n_modes,
        ref grid,
        y_range,
        port_modes,
        port_z0_ohm,
        reconstruct_fields,
    } = problem.study
    else {
        return Err("assemble_eme requires Study::Eme".into());
    };
    if !(freq > 0.0 && freq.is_finite()) {
        return Err(format!("EME frequency must be positive and finite (got {freq})"));
    }
    if n_modes == 0 {
        return Err("EME requires n_modes >= 1".into());
    }
    validate_grid_spec(grid)?;

    let d3 = discretise(mesh).map_err(|e| e.to_string())?;
    let (lo3, hi3) = d3.bbox();
    let [y0, y1] = y_range.unwrap_or([lo3[1], hi3[1]]);
    let span_y = y1 - y0;
    if !(span_y > 0.0 && span_y.is_finite()) {
        return Err(format!("EME propagation interval [y0, y1] = [{y0}, {y1}] must satisfy y1 > y0"));
    }
    let tol_y = 1e-7 * span_y.max(1e-12);
    if y0 < lo3[1] - tol_y || y1 > hi3[1] + tol_y {
        return Err(format!(
            "EME y_range [{y0:.6e}, {y1:.6e}] lies outside the 3D mesh y-bounds [{:.6e}, {:.6e}]",
            lo3[1], hi3[1]
        ));
    }

    // Reject PML absorption along propagation axis y on cross-section slices
    if let Some((dom, _)) = problem.materials.pml.iter().find(|(_, p)| p.y != femcore::props::PmlDir::Off) {
        return Err(format!(
            "EME domain {dom} has a PML absorbing along propagation axis y; EME handles propagation along +y analytically and uses transverse (x, z) PMLs only"
        ));
    }

    // Build segment plans and cascade tree
    let mut raw_segments: Vec<EmeSegmentPlan> = Vec::new();
    let mut subgrid_counter = 0usize;
    let tree = build_grid_tree(grid, y0, y1, hi3[1], n_modes, &mut raw_segments, &mut subgrid_counter)?;
    if raw_segments.is_empty() {
        return Err("EME grid produced zero segments".into());
    }

    // Prepare 3D temperature and potential transfer evaluators (cross-mesh FR-44 or same-mesh)
    let temp_src_disc_owned: Option<std::sync::Arc<Discretisation>> =
        if let Some(ref tt) = problem.transferred_temperature {
            match &tt.0.discretisation {
                Some(d) => Some(d.clone()),
                None => Some(std::sync::Arc::new(discretise(&tt.0.mesh).map_err(|e| e.to_string())?)),
            }
        } else {
            None
        };
    let temp_idx_and_field: Option<(TransferIndex<'_>, SourceField<'_>)> = if let (Some(ref tt), Some(ref sd)) =
        (&problem.transferred_temperature, &temp_src_disc_owned)
    {
        let idx = TransferIndex::new(SourceMesh::Simplex(sd.as_ref()), &tt.0.domain_names, TransferOptions::default())
            .map_err(|e| e.to_string())?;
        let sf =
            SourceField::Nodal { values: &tt.0.values, ncomp: tt.0.components.max(1), order: tt.0.order.clamp(1, 2) };
        Some((idx, sf))
    } else if let Some(ref tf_vec) = problem.temperature_field {
        if !tf_vec.is_empty() {
            let idx = TransferIndex::new(SourceMesh::Simplex(&d3), &problem.domain_names, TransferOptions::default())
                .map_err(|e| e.to_string())?;
            let sf = if tf_vec.len() == d3.n_cells() && tf_vec.len() != d3.n_points() {
                SourceField::PerCell { values: tf_vec }
            } else if tf_vec.len() == d3.n_dofs(2) {
                SourceField::Nodal { values: tf_vec, ncomp: 1, order: 2 }
            } else {
                SourceField::Nodal { values: tf_vec, ncomp: 1, order: 1 }
            };
            Some((idx, sf))
        } else {
            None
        }
    } else if let Some(ref ct_vec) = problem.cell_temperatures {
        if ct_vec.len() == d3.n_cells() {
            let idx = TransferIndex::new(SourceMesh::Simplex(&d3), &problem.domain_names, TransferOptions::default())
                .map_err(|e| e.to_string())?;
            Some((idx, SourceField::PerCell { values: ct_vec }))
        } else {
            None
        }
    } else {
        None
    };

    let pot_src_disc_owned: Option<std::sync::Arc<Discretisation>> = if let Some(ref tp) = problem.transferred_potential
    {
        match &tp.0.discretisation {
            Some(d) => Some(d.clone()),
            None => Some(std::sync::Arc::new(discretise(&tp.0.mesh).map_err(|e| e.to_string())?)),
        }
    } else {
        None
    };
    let pot_idx_and_field: Option<(TransferIndex<'_>, SourceField<'_>)> = if !problem.sheet_2d.is_empty() {
        if let (Some(ref tp), Some(ref sd)) = (&problem.transferred_potential, &pot_src_disc_owned) {
            let idx =
                TransferIndex::new(SourceMesh::Simplex(sd.as_ref()), &tp.0.domain_names, TransferOptions::default())
                    .map_err(|e| e.to_string())?;
            let sf = SourceField::Nodal {
                values: &tp.0.values,
                ncomp: tp.0.components.max(1),
                order: tp.0.order.clamp(1, 2),
            };
            Some((idx, sf))
        } else if let Some(ref pv) = problem.imported_potential {
            if !pv.is_empty() {
                let idx =
                    TransferIndex::new(SourceMesh::Simplex(&d3), &problem.domain_names, TransferOptions::default())
                        .map_err(|e| e.to_string())?;
                let ord = if pv.len() == d3.n_dofs(2) { 2 } else { 1 };
                Some((idx, SourceField::Nodal { values: pv, ncomp: 1, order: ord }))
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    };

    // Precompute sorted unique y-coordinates of the 3D mesh for fast slab detection and vertex avoidance
    let mesh_y_planes = collect_sorted_y_planes(&d3, span_y);
    let quad4 = quad_rule(2, 4);
    let quad6 = quad_rule(2, 6);

    // Prepare 3D/2D electro-optic DC electric-field and photoelastic strain transfer evaluators (FR-37 + FR-39)
    let eo_pe_domains = problem.domain_eo_pe_global();
    let has_eo_dom = eo_pe_domains.values().any(|(r, _)| r.is_some());
    let has_pe_dom = eo_pe_domains.values().any(|(_, p)| p.is_some());

    let elec_src_disc_owned: Option<std::sync::Arc<Discretisation>> = if has_eo_dom {
        if let Some(ref te) = problem.transferred_electric_field {
            match &te.0.discretisation {
                Some(d) => Some(d.clone()),
                None => Some(std::sync::Arc::new(discretise(&te.0.mesh).map_err(|e| e.to_string())?)),
            }
        } else {
            None
        }
    } else {
        None
    };
    let mut elec_owned_cells: Vec<[f64; 3]> = Vec::new();
    if has_eo_dom && elec_src_disc_owned.is_none() {
        if let Some(ref ef_vec) = problem.electric_field {
            if ef_vec.len() == 3 * d3.n_cells() && ef_vec.len() != 3 * d3.n_points() {
                elec_owned_cells =
                    (0..d3.n_cells()).map(|c| [ef_vec[3 * c], ef_vec[3 * c + 1], ef_vec[3 * c + 2]]).collect();
            }
        }
    }
    let elec_idx_and_src: Option<(TransferIndex<'_>, VecSliceSource<'_, 3>)> = if has_eo_dom {
        if let (Some(ref te), Some(ref sd)) = (&problem.transferred_electric_field, &elec_src_disc_owned) {
            let idx =
                TransferIndex::new(SourceMesh::Simplex(sd.as_ref()), &te.0.domain_names, TransferOptions::default())
                    .map_err(|e| e.to_string())?;
            let ncomp = te.0.components.max(1);
            let sf = SourceField::Nodal { values: &te.0.values, ncomp, order: te.0.order.clamp(1, 2) };
            Some((idx, VecSliceSource::Nodal { sf, src_dim: sd.dim, ncomp }))
        } else if !elec_owned_cells.is_empty() {
            let idx = TransferIndex::new(SourceMesh::Simplex(&d3), &problem.domain_names, TransferOptions::default())
                .map_err(|e| e.to_string())?;
            Some((idx, VecSliceSource::PerCell3d(&elec_owned_cells)))
        } else if let Some(ce_vec) = problem.cell_electric_field.as_ref().filter(|v| v.len() == d3.n_cells()) {
            let idx = TransferIndex::new(SourceMesh::Simplex(&d3), &problem.domain_names, TransferOptions::default())
                .map_err(|e| e.to_string())?;
            Some((idx, VecSliceSource::PerCell3d(ce_vec)))
        } else if let Some(ef_vec) = problem.electric_field.as_ref().filter(|v| !v.is_empty()) {
            let ord = if ef_vec.len() == 3 * d3.n_dofs(2) {
                2
            } else if ef_vec.len() == 3 * d3.n_points() {
                1
            } else {
                0
            };
            if ord > 0 {
                let idx =
                    TransferIndex::new(SourceMesh::Simplex(&d3), &problem.domain_names, TransferOptions::default())
                        .map_err(|e| e.to_string())?;
                Some((
                    idx,
                    VecSliceSource::Nodal {
                        sf: SourceField::Nodal { values: ef_vec, ncomp: 3, order: ord },
                        src_dim: 3,
                        ncomp: 3,
                    },
                ))
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    };

    let strain_src_disc_owned: Option<std::sync::Arc<Discretisation>> = if has_pe_dom {
        if let Some(ref ts) = problem.transferred_strain {
            match &ts.0.discretisation {
                Some(d) => Some(d.clone()),
                None => Some(std::sync::Arc::new(discretise(&ts.0.mesh).map_err(|e| e.to_string())?)),
            }
        } else {
            None
        }
    } else {
        None
    };
    let mut strain_owned_cells: Vec<[f64; 6]> = Vec::new();
    if has_pe_dom && strain_src_disc_owned.is_none() {
        if let Some(ref sf_vec) = problem.strain_field {
            if sf_vec.len() == 6 * d3.n_cells() && sf_vec.len() != 6 * d3.n_points() {
                strain_owned_cells = (0..d3.n_cells())
                    .map(|c| {
                        [
                            sf_vec[6 * c],
                            sf_vec[6 * c + 1],
                            sf_vec[6 * c + 2],
                            sf_vec[6 * c + 3],
                            sf_vec[6 * c + 4],
                            sf_vec[6 * c + 5],
                        ]
                    })
                    .collect();
            }
        }
    }
    let strain_idx_and_src: Option<(TransferIndex<'_>, VecSliceSource<'_, 6>)> = if has_pe_dom {
        if let (Some(ref ts), Some(ref sd)) = (&problem.transferred_strain, &strain_src_disc_owned) {
            let idx =
                TransferIndex::new(SourceMesh::Simplex(sd.as_ref()), &ts.0.domain_names, TransferOptions::default())
                    .map_err(|e| e.to_string())?;
            let ncomp = ts.0.components.max(1);
            let sf = SourceField::Nodal { values: &ts.0.values, ncomp, order: ts.0.order.clamp(1, 2) };
            Some((idx, VecSliceSource::Nodal { sf, src_dim: sd.dim, ncomp }))
        } else if !strain_owned_cells.is_empty() {
            let idx = TransferIndex::new(SourceMesh::Simplex(&d3), &problem.domain_names, TransferOptions::default())
                .map_err(|e| e.to_string())?;
            Some((idx, VecSliceSource::PerCell3d(&strain_owned_cells)))
        } else if let Some(cs_vec) = problem.cell_strain.as_ref().filter(|v| v.len() == d3.n_cells()) {
            let idx = TransferIndex::new(SourceMesh::Simplex(&d3), &problem.domain_names, TransferOptions::default())
                .map_err(|e| e.to_string())?;
            Some((idx, VecSliceSource::PerCell3d(cs_vec)))
        } else if let Some(sf_vec) = problem.strain_field.as_ref().filter(|v| !v.is_empty()) {
            let ord = if sf_vec.len() == 6 * d3.n_dofs(2) {
                2
            } else if sf_vec.len() == 6 * d3.n_points() {
                1
            } else {
                0
            };
            if ord > 0 {
                let idx =
                    TransferIndex::new(SourceMesh::Simplex(&d3), &problem.domain_names, TransferOptions::default())
                        .map_err(|e| e.to_string())?;
                Some((
                    idx,
                    VecSliceSource::Nodal {
                        sf: SourceField::Nodal { values: sf_vec, ncomp: 6, order: ord },
                        src_dim: 3,
                        ncomp: 6,
                    },
                ))
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    };

    let mut unique_slices: Vec<EmeUniqueSlice> = Vec::new();
    let mut by_hash: HashMap<String, Vec<usize>> = HashMap::new();

    for seg in &mut raw_segments {
        femcore::progress::check()?;
        let (mesh_2d, y_actual) = extract_2d_cross_section(&d3, &mesh_y_planes, seg.y_mesh_lo, seg.y_mesh_hi)?;
        seg.y_sample = y_actual;

        if let Some(r_bend) = seg.bend_radius_m {
            for v in 0..mesh_2d.n_points() {
                let x = mesh_2d.point(v)[0];
                let _ = conformal_bend_coord(x, r_bend)?;
            }
        }

        // Sample domain-filtered temperatures T(x_q, y_c, z_q) when thermo-optic coupling is active
        let (cell_quad4_temps, cell_quad6_temps, cell_avg_temps, has_subcell_temp_grad) =
            if !problem.thermo_optic_domains.is_empty() && temp_idx_and_field.is_some() {
                let (ref t_idx, ref t_sf) =
                    temp_idx_and_field.as_ref().ok_or_else(|| "internal temperature index missing".to_string())?;
                sample_slice_temperatures(&mesh_2d, seg.y_sample, problem, t_idx, t_sf, &quad4, &quad6)?
            } else {
                (Vec::new(), Vec::new(), Vec::new(), false)
            };

        let (cell_quad4_e_dc, cell_quad6_e_dc, cell_avg_e_dc, subcell_e_grad) =
            if let Some((ref e_idx, ref e_src)) = elec_idx_and_src {
                sample_slice_vec_field::<3>(&mesh_2d, seg.y_sample, problem, e_idx, e_src, &quad4, &quad6)?
            } else {
                (Vec::new(), Vec::new(), Vec::new(), false)
            };

        let (cell_quad4_strain, cell_quad6_strain, cell_avg_strain, subcell_s_grad) =
            if let Some((ref s_idx, ref s_src)) = strain_idx_and_src {
                sample_slice_vec_field::<6>(&mesh_2d, seg.y_sample, problem, s_idx, s_src, &quad4, &quad6)?
            } else {
                (Vec::new(), Vec::new(), Vec::new(), false)
            };
        let has_subcell_eo_pe_grad = subcell_e_grad || subcell_s_grad;

        let slice_potential = if let Some((ref p_idx, ref p_sf)) = pot_idx_and_field {
            Some(sample_slice_potential(&mesh_2d, seg.y_sample, p_idx, p_sf)?)
        } else {
            None
        };

        let hash = compute_slice_blake3(
            &mesh_2d,
            problem,
            freq,
            seg.bend_radius_m,
            &eo_pe_domains,
            &cell_quad4_temps,
            &cell_quad4_e_dc,
            &cell_quad4_strain,
            slice_potential.as_deref(),
        )?;

        let matched_uid = by_hash.get(&hash).and_then(|cands| {
            cands.iter().copied().find(|&uid| slices_geometrically_compatible(&unique_slices[uid].mesh_2d, &mesh_2d))
        });

        if let Some(uid) = matched_uid {
            seg.unique_slice = uid;
            if seg.n_modes > unique_slices[uid].max_modes {
                unique_slices[uid].max_modes = seg.n_modes;
            }
        } else {
            let uid = unique_slices.len();
            seg.unique_slice = uid;
            by_hash.entry(hash.clone()).or_default().push(uid);
            unique_slices.push(EmeUniqueSlice {
                unique_id: uid,
                blake3_hash: hash,
                sample_y_m: seg.y_sample,
                mesh_2d,
                max_modes: seg.n_modes,
                bend_radius_m: seg.bend_radius_m,
                cell_quad4_temps,
                cell_quad6_temps,
                cell_avg_temps,
                has_subcell_temp_grad,
                cell_quad4_e_dc,
                cell_quad6_e_dc,
                cell_avg_e_dc,
                cell_quad4_strain,
                cell_quad6_strain,
                cell_avg_strain,
                has_subcell_eo_pe_grad,
                slice_potential,
            });
        }
    }

    // Update segment unique_slice references inside EmeCascadeNode tree
    let mut tree = tree;
    assign_tree_unique_slices(&mut tree, &raw_segments);

    // Extract optional explicit waveports on y = y0 and y = y1
    let waveport_patches = extract_eme_waveport_patches(&d3, &problem.waveports, &problem.pmc, y0, y1)?;

    Ok(EmeAssembly {
        mesh_3d: d3,
        problem: problem.clone(),
        freq_hz: freq,
        n_target,
        evanescent,
        y_range_m: [y0, y1],
        port_modes,
        port_z0_ohm,
        reconstruct_fields,
        tree,
        total_periodic_subgrids: subgrid_counter,
        mesh_segments: raw_segments,
        unique_slices,
        waveport_patches,
    })
}

fn build_grid_tree(
    spec: &EmeGridSpec,
    y_lo: f64,
    y_hi: f64,
    mesh_y_max: f64,
    default_modes: usize,
    segments: &mut Vec<EmeSegmentPlan>,
    subgrid_counter: &mut usize,
) -> Result<EmeCascadeNode, String> {
    let span = y_hi - y_lo;
    if !(span > 0.0 && span.is_finite()) {
        return Err(format!("EME subgrid interval [{y_lo:.6e}, {y_hi:.6e}] has non-positive span"));
    }
    match spec {
        EmeGridSpec::Uniform { n_cells, n_modes, bend_radius } => {
            let nc = (*n_cells).max(1);
            let m = n_modes.unwrap_or(default_modes).max(1);
            let dl = span / (nc as f64);
            let mut children = Vec::with_capacity(nc);
            for k in 0..nc {
                let s_idx = segments.len();
                let a = y_lo + (k as f64) * dl;
                let b = if k + 1 == nc { y_hi } else { y_lo + ((k + 1) as f64) * dl };
                segments.push(EmeSegmentPlan {
                    segment_index: s_idx,
                    y_mesh_lo: a,
                    y_mesh_hi: b,
                    y_sample: 0.5 * (a + b),
                    length_m: dl,
                    n_modes: m,
                    bend_radius_m: *bend_radius,
                    unique_slice: 0,
                });
                children.push(EmeCascadeNode::Segment {
                    segment_index: s_idx,
                    unique_slice: 0,
                    length_m: dl,
                    n_modes: m,
                    bend_radius_m: *bend_radius,
                });
            }
            if children.len() == 1 {
                Ok(children.remove(0))
            } else {
                Ok(EmeCascadeNode::Composite { children })
            }
        }
        EmeGridSpec::Explicit { lengths, n_modes, bend_radius, mesh_lengths } => {
            let l_sum: f64 = lengths.iter().sum();
            if !(l_sum > 0.0 && l_sum.is_finite()) {
                return Err("EME Explicit subgrid total length must be positive".into());
            }
            let m_lens = mesh_lengths.as_deref().unwrap_or(lengths);
            let ml_sum: f64 = m_lens.iter().sum();
            if !(ml_sum > 0.0 && ml_sum.is_finite()) {
                return Err("EME Explicit subgrid mesh_lengths sum must be positive".into());
            }
            // Map mesh_lengths (or lengths) proportionally onto mesh interval [y_lo, y_hi] while keeping true physical `length_m = lengths[i]`
            let scale = span / ml_sum;
            let mut cursor = y_lo;
            let mut children = Vec::with_capacity(lengths.len());
            for (i, &len) in lengths.iter().enumerate() {
                let s_idx = segments.len();
                let dy_mesh = m_lens[i] * scale;
                let a = cursor;
                let b = if i + 1 == lengths.len() { y_hi } else { cursor + dy_mesh };
                cursor = b;
                let m = n_modes.as_ref().map_or(default_modes, |ms| ms[i]).max(1);
                let br = bend_radius.as_ref().and_then(|bs| bs[i]);
                segments.push(EmeSegmentPlan {
                    segment_index: s_idx,
                    y_mesh_lo: a,
                    y_mesh_hi: b,
                    y_sample: 0.5 * (a + b),
                    length_m: len,
                    n_modes: m,
                    bend_radius_m: br,
                    unique_slice: 0,
                });
                children.push(EmeCascadeNode::Segment {
                    segment_index: s_idx,
                    unique_slice: 0,
                    length_m: len,
                    n_modes: m,
                    bend_radius_m: br,
                });
            }
            if children.len() == 1 {
                Ok(children.remove(0))
            } else {
                Ok(EmeCascadeNode::Composite { children })
            }
        }
        EmeGridSpec::PeriodicSubgrid { unit_cell, num_reps, unit_length } => {
            let sg_idx = *subgrid_counter;
            *subgrid_counter += 1;
            let reps = (*num_reps).max(1);
            // Determine unit-cell mesh sampling interval [y_lo, y_uc_hi]
            let uc_opt_len = unit_length.or_else(|| unit_cell.explicit_unit_mesh_length());
            let (uc_mesh_hi, uc_override_scale) = if let Some(lu) = uc_opt_len {
                let tol = 1e-6 * span.max(lu);
                if (span - lu).abs() <= tol || y_lo + lu <= mesh_y_max + tol {
                    (y_lo + lu.min(span), None)
                } else {
                    (y_hi, Some(lu))
                }
            } else {
                // Unit cell is Uniform without explicit length: if mesh covers `reps` periods, sample first period `span / reps`
                (y_lo + span / (reps as f64), None)
            };
            let seg_start = segments.len();
            let mut uc_node =
                build_grid_tree(unit_cell, y_lo, uc_mesh_hi, mesh_y_max, default_modes, segments, subgrid_counter)?;
            if let Some(lu) = uc_override_scale {
                let cur_sum: f64 = segments[seg_start..].iter().map(|s| s.length_m).sum();
                if cur_sum > 0.0 {
                    let ratio = lu / cur_sum;
                    for s in &mut segments[seg_start..] {
                        s.length_m *= ratio;
                    }
                    scale_tree_lengths(&mut uc_node, ratio);
                }
            } else if let Some(lu) = *unit_length {
                // Ensure Uniform unit_cell with explicit unit_length propagates over `lu`
                let cur_sum: f64 = segments[seg_start..].iter().map(|s| s.length_m).sum();
                if cur_sum > 0.0 && (cur_sum - lu).abs() > 1e-15 * lu {
                    let ratio = lu / cur_sum;
                    for s in &mut segments[seg_start..] {
                        s.length_m *= ratio;
                    }
                    scale_tree_lengths(&mut uc_node, ratio);
                }
            }
            Ok(EmeCascadeNode::Periodic { subgrid_index: sg_idx, unit_cell: Box::new(uc_node), num_reps: reps })
        }
        EmeGridSpec::Composite { subgrids } => {
            // Allocate mesh y-intervals to child subgrids: check full optical length vs compact 1-period mesh length
            let full_exp: f64 = subgrids.iter().filter_map(EmeGridSpec::explicit_optical_length).sum();
            let compact_exp: f64 = subgrids.iter().filter_map(EmeGridSpec::explicit_unit_mesh_length).sum();
            let uniform_weight: usize = subgrids.iter().map(EmeGridSpec::uniform_cell_weight).sum();
            let tol = 1e-5 * span.max(1e-12);
            let use_compact_periodic = full_exp > span + tol && compact_exp <= span + tol;
            let fixed_sum = if use_compact_periodic { compact_exp } else { full_exp };
            let rem_span = (span - fixed_sum).max(0.0);

            let mut cursor = y_lo;
            let mut children = Vec::with_capacity(subgrids.len());
            for (i, sg) in subgrids.iter().enumerate() {
                let sg_span = if let Some(len) =
                    if use_compact_periodic { sg.explicit_unit_mesh_length() } else { sg.explicit_optical_length() }
                {
                    if uniform_weight == 0 && fixed_sum > 0.0 {
                        len * (span / fixed_sum)
                    } else {
                        len
                    }
                } else {
                    let w = sg.uniform_cell_weight().max(1);
                    rem_span * (w as f64) / (uniform_weight.max(1) as f64)
                };
                let b = if i + 1 == subgrids.len() { y_hi } else { (cursor + sg_span).min(y_hi) };
                let node = build_grid_tree(sg, cursor, b, mesh_y_max, default_modes, segments, subgrid_counter)?;
                cursor = b;
                children.push(node);
            }
            Ok(EmeCascadeNode::Composite { children })
        }
    }
}

fn scale_tree_lengths(node: &mut EmeCascadeNode, ratio: f64) {
    match node {
        EmeCascadeNode::Segment { length_m, .. } => *length_m *= ratio,
        EmeCascadeNode::Composite { children } => {
            for ch in children {
                scale_tree_lengths(ch, ratio);
            }
        }
        EmeCascadeNode::Periodic { unit_cell, .. } => scale_tree_lengths(unit_cell, ratio),
    }
}

fn assign_tree_unique_slices(node: &mut EmeCascadeNode, segments: &[EmeSegmentPlan]) {
    match node {
        EmeCascadeNode::Segment { segment_index, unique_slice, .. } => {
            *unique_slice = segments[*segment_index].unique_slice;
        }
        EmeCascadeNode::Composite { children } => {
            for ch in children {
                assign_tree_unique_slices(ch, segments);
            }
        }
        EmeCascadeNode::Periodic { unit_cell, .. } => assign_tree_unique_slices(unit_cell, segments),
    }
}

fn collect_sorted_y_planes(d3: &SimplexMesh, span_y: f64) -> Vec<f64> {
    let mut ys: Vec<f64> = (0..d3.n_points()).map(|v| d3.point(v)[1]).collect();
    ys.sort_by(f64::total_cmp);
    let tol = 1e-10 * span_y.max(1e-12);
    let mut planes: Vec<f64> = Vec::new();
    for y in ys {
        if planes.last().is_none_or(|&prev| (y - prev).abs() > tol) {
            planes.push(y);
        }
    }
    planes
}

/// Extracts a conforming 2D `(x, z)` `SimplexMesh` from `d3` inside segment `[y_lo, y_hi]`.
///
/// 1. **Exact vertical-prism slab fast path**: when the 3D tetrahedra straddling `y_mid` form a straight
///    vertical slab `[Y_k, Y_{k+1}]` (as produced by `extrude_y`, `sweep_model_3d`, or `block`), extracts
///    the exact un-subdivided 2D base triangulation on `Y_k` — zero extra Steiner vertices from diagonal
///    tet cuts, identical 2D mesh at every `y` inside uniform waveguides, and machine-precision agreement
///    with 2D `ems::modes`.
/// 2. **General marching-tetrahedra path**: for sloped sidewalls (adiabatic tapers) or unstructured 3D meshes,
///    chooses `y_sample in [y_lo + 0.15 L, y_hi - 0.15 L]` strictly away from all 3D mesh vertices and
///    intersects every straddling tetrahedron and 3D boundary/sheet facet along sorted 3D edges.
fn extract_2d_cross_section(
    d3: &SimplexMesh,
    y_planes: &[f64],
    y_lo: f64,
    y_hi: f64,
) -> Result<(SimplexMesh, f64), String> {
    let y_mid = 0.5 * (y_lo + y_hi);
    let seg_len = (y_hi - y_lo).max(1e-15);

    // Choose a vertex-free plane `y_cut` near `y_mid` inside `[y_lo + 0.15 L, y_hi - 0.15 L]`
    let win_lo = y_lo + 0.15 * seg_len;
    let win_hi = y_hi - 0.15 * seg_len;
    let min_clear = 1e-4 * seg_len;
    let near_mid = y_planes.iter().any(|&yp| (yp - y_mid).abs() <= min_clear);
    let y_cut = if !near_mid {
        y_mid
    } else {
        let mut pts_in = vec![win_lo];
        for &yp in y_planes {
            if yp > win_lo && yp < win_hi {
                pts_in.push(yp);
            }
        }
        pts_in.push(win_hi);
        let mut best_y = y_mid;
        let mut best_gap = -1.0_f64;
        for w in pts_in.windows(2) {
            let gap = w[1] - w[0];
            let cand = 0.5 * (w[0] + w[1]);
            if gap > best_gap + 1e-15 * seg_len
                || ((gap - best_gap).abs() <= 1e-15 * seg_len && (cand - y_mid).abs() < (best_y - y_mid).abs())
            {
                best_gap = gap;
                best_y = cand;
            }
        }
        best_y
    };

    if let Some(slab_mesh) = try_extract_vertical_prism_slab(d3, y_planes, y_cut)? {
        return Ok((canonicalize_2d_mesh(&slab_mesh), y_mid));
    }

    let sliced = slice_tetrahedra_at_plane(d3, y_cut)?;
    Ok((canonicalize_2d_mesh(&sliced), y_cut))
}

/// Checks whether the tetrahedra straddling `y_cut` form a straight vertical prism slab `[y_a, y_b]`
/// and, if so, extracts the exact 2D base triangulation on `y_a`.
fn try_extract_vertical_prism_slab(
    d3: &SimplexMesh,
    y_planes: &[f64],
    y_cut: f64,
) -> Result<Option<SimplexMesh>, String> {
    if y_planes.len() < 2 {
        return Ok(None);
    }
    let idx_hi = match y_planes.binary_search_by(|yp| yp.total_cmp(&y_cut)) {
        Ok(_) => return Ok(None),
        Err(i) => i,
    };
    if idx_hi == 0 || idx_hi >= y_planes.len() {
        return Ok(None);
    }
    let y_a = y_planes[idx_hi - 1];
    let y_b = y_planes[idx_hi];
    let (lo3, hi3) = d3.bbox();
    let diag_xz = (hi3[0] - lo3[0]).hypot(hi3[2] - lo3[2]).max(1e-12);
    let tol_y = 1e-9 * (hi3[1] - lo3[1]).max(1e-12);
    let tol_xz = 1e-9 * diag_xz;

    let mut straddling_cells = Vec::new();
    for ci in 0..d3.n_cells() {
        let cv = d3.cell(ci);
        let mut min_y = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        for &v in cv {
            let yv = d3.point(v as usize)[1];
            min_y = min_y.min(yv);
            max_y = max_y.max(yv);
        }
        if min_y < y_cut && max_y > y_cut {
            // Every vertex of every straddling cell must sit on y_a or y_b
            for &v in cv {
                let yv = d3.point(v as usize)[1];
                if (yv - y_a).abs() > tol_y && (yv - y_b).abs() > tol_y {
                    return Ok(None);
                }
            }
            straddling_cells.push(ci);
        }
    }
    if straddling_cells.is_empty() {
        return Ok(None);
    }

    // Map quantized (x, z) on y_a -> 2D vertex index, and verify every y_b vertex in straddling cells has an exact vertical partner on y_a!
    let q_scale = 1e10 / diag_xz;
    let q_key = |x: f64, z: f64| -> (i64, i64) {
        (((x - lo3[0]) * q_scale).round() as i64, ((z - lo3[2]) * q_scale).round() as i64)
    };

    let mut ya_map: HashMap<(i64, i64), u32> = HashMap::new();
    let mut v3_to_v2: HashMap<u32, u32> = HashMap::new();
    let mut pts_2d: Vec<f64> = Vec::new();

    // First register all vertices on y_a touched by straddling cells
    for &ci in &straddling_cells {
        for &v in d3.cell(ci) {
            let p = d3.point(v as usize);
            if (p[1] - y_a).abs() <= tol_y {
                if let std::collections::hash_map::Entry::Vacant(e) = v3_to_v2.entry(v) {
                    let key = q_key(p[0], p[2]);
                    let v2 = *ya_map.entry(key).or_insert_with(|| {
                        let id = (pts_2d.len() / 2) as u32;
                        pts_2d.push(p[0]);
                        pts_2d.push(p[2]);
                        id
                    });
                    e.insert(v2);
                }
            }
        }
    }

    // Verify every y_b vertex in straddling cells projects vertically onto an existing y_a vertex within tol_xz
    for &ci in &straddling_cells {
        for &v in d3.cell(ci) {
            let p = d3.point(v as usize);
            if (p[1] - y_b).abs() <= tol_y && !v3_to_v2.contains_key(&v) {
                let key = q_key(p[0], p[2]);
                let mut matched = ya_map.get(&key).copied();
                if matched.is_none() {
                    for dx in -1..=1 {
                        for dz in -1..=1 {
                            if let Some(&cand) = ya_map.get(&(key.0 + dx, key.1 + dz)) {
                                matched = Some(cand);
                                break;
                            }
                        }
                    }
                }
                let Some(v2) = matched else {
                    return Ok(None);
                };
                let px = pts_2d[2 * v2 as usize];
                let pz = pts_2d[2 * v2 as usize + 1];
                if (p[0] - px).hypot(p[2] - pz) > tol_xz {
                    return Ok(None);
                }
                v3_to_v2.insert(v, v2);
            }
        }
    }

    // Collect base triangles (cells with 3 vertices on y_a and 1 on y_b) and check consistent prism labels
    let mut tri_label_map: BTreeMap<[u32; 3], (u32, u32, u32, i32)> = BTreeMap::new();
    for &ci in &straddling_cells {
        let cv = d3.cell(ci);
        let mut bot = Vec::with_capacity(3);
        for &v in cv {
            if (d3.point(v as usize)[1] - y_a).abs() <= tol_y {
                bot.push(v3_to_v2[&v]);
            }
        }
        // Every straddling tet of a vertical prism projects onto the prism's 3 unique base vertices!
        let mut proj: Vec<u32> = cv.iter().map(|v| v3_to_v2[v]).collect();
        proj.sort_unstable();
        proj.dedup();
        if proj.len() != 3 {
            return Ok(None);
        }
        let sorted_key = [proj[0], proj[1], proj[2]];
        let lab = d3.cell_labels[ci];
        if let Some(existing) = tri_label_map.get(&sorted_key) {
            if existing.3 != lab {
                return Ok(None);
            }
        } else if bot.len() == 3 {
            tri_label_map.insert(sorted_key, (bot[0], bot[1], bot[2], lab));
        }
    }
    if tri_label_map.is_empty() {
        return Ok(None);
    }

    let mut cells_2d = Vec::with_capacity(3 * tri_label_map.len());
    let mut cell_labels_2d = Vec::with_capacity(tri_label_map.len());
    for (_, (mut a, mut b, mut c_v, lab)) in tri_label_map {
        orient_ccw_2d(&pts_2d, &mut a, &mut b, &mut c_v);
        cells_2d.extend_from_slice(&[a, b, c_v]);
        cell_labels_2d.push(lab);
    }

    // Extract 1D facets from 3D facets spanning [y_a, y_b] (2 vertices on y_a, 1 vertex on y_b)
    let mut edge_label_map: BTreeMap<[u32; 2], i32> = BTreeMap::new();
    for f in 0..d3.n_facets() {
        let fv = d3.facet(f);
        let mut on_a = Vec::with_capacity(3);
        let mut on_b = 0usize;
        for &v in fv {
            let yv = d3.point(v as usize)[1];
            if (yv - y_a).abs() <= tol_y {
                if let Some(&v2) = v3_to_v2.get(&v) {
                    on_a.push(v2);
                }
            } else if (yv - y_b).abs() <= tol_y && v3_to_v2.contains_key(&v) {
                on_b += 1;
            }
        }
        if on_a.len() == 2 && on_b == 1 && on_a[0] != on_a[1] {
            let key = [on_a[0].min(on_a[1]), on_a[0].max(on_a[1])];
            edge_label_map.entry(key).or_insert(d3.facet_labels[f]);
        }
    }

    let mut facets_2d = Vec::with_capacity(2 * edge_label_map.len());
    let mut facet_labels_2d = Vec::with_capacity(edge_label_map.len());
    for ([u, w], lab) in edge_label_map {
        facets_2d.extend_from_slice(&[u, w]);
        facet_labels_2d.push(lab);
    }

    let mesh = SimplexMesh::new(2, pts_2d, cells_2d, cell_labels_2d, facets_2d, facet_labels_2d);
    if mesh.check().is_err() {
        return Ok(None);
    }
    Ok(Some(mesh))
}

fn edge_cut_vertex(
    d3: &SimplexMesh,
    y_cut: f64,
    edge_to_v2: &mut HashMap<(u32, u32), u32>,
    pts_2d: &mut Vec<f64>,
    u: u32,
    w: u32,
) -> u32 {
    let key = (u.min(w), u.max(w));
    *edge_to_v2.entry(key).or_insert_with(|| {
        let pu = d3.point(key.0 as usize);
        let pw = d3.point(key.1 as usize);
        let denom = pw[1] - pu[1];
        let t = if denom.abs() > 1e-300 { ((y_cut - pu[1]) / denom).clamp(0.0, 1.0) } else { 0.5 };
        let x = pu[0] + t * (pw[0] - pu[0]);
        let z = pu[2] + t * (pw[2] - pu[2]);
        let id = (pts_2d.len() / 2) as u32;
        pts_2d.push(x);
        pts_2d.push(z);
        id
    })
}

/// General edge-based marching-tetrahedra slice at `y = y_cut` (`y_cut` strictly away from all 3D vertices).
fn slice_tetrahedra_at_plane(d3: &SimplexMesh, y_cut: f64) -> Result<SimplexMesh, String> {
    let mut edge_to_v2: HashMap<(u32, u32), u32> = HashMap::new();
    let mut pts_2d: Vec<f64> = Vec::new();

    let mut cells_2d: Vec<u32> = Vec::new();
    let mut cell_labels_2d: Vec<i32> = Vec::new();

    for ci in 0..d3.n_cells() {
        let cv = d3.cell(ci);
        let mut neg = Vec::with_capacity(4);
        let mut pos = Vec::with_capacity(4);
        for &v in cv {
            if d3.point(v as usize)[1] < y_cut {
                neg.push(v);
            } else {
                pos.push(v);
            }
        }
        if neg.is_empty() || pos.is_empty() {
            continue;
        }
        let lab = d3.cell_labels[ci];
        if neg.len() == 1 && pos.len() == 3 {
            let u = neg[0];
            let (mut a, mut b, mut c_v) = (
                edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, u, pos[0]),
                edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, u, pos[1]),
                edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, u, pos[2]),
            );
            orient_ccw_2d(&pts_2d, &mut a, &mut b, &mut c_v);
            cells_2d.extend_from_slice(&[a, b, c_v]);
            cell_labels_2d.push(lab);
        } else if neg.len() == 3 && pos.len() == 1 {
            let w = pos[0];
            let (mut a, mut b, mut c_v) = (
                edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[0], w),
                edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[1], w),
                edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[2], w),
            );
            orient_ccw_2d(&pts_2d, &mut a, &mut b, &mut c_v);
            cells_2d.extend_from_slice(&[a, b, c_v]);
            cell_labels_2d.push(lab);
        } else {
            // (2 vs 2) convex quad: cyclic perimeter is p00 -> p10 -> p11 -> p01
            neg.sort_unstable();
            pos.sort_unstable();
            let p00 = edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[0], pos[0]);
            let p10 = edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[1], pos[0]);
            let p11 = edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[1], pos[1]);
            let p01 = edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[0], pos[1]);

            let d0_sq = dist2_2d(&pts_2d, p00, p11);
            let d1_sq = dist2_2d(&pts_2d, p10, p01);
            let tris =
                if d0_sq <= d1_sq { [(p00, p10, p11), (p00, p11, p01)] } else { [(p00, p10, p01), (p10, p11, p01)] };
            for (mut a, mut b, mut c_v) in tris {
                orient_ccw_2d(&pts_2d, &mut a, &mut b, &mut c_v);
                cells_2d.extend_from_slice(&[a, b, c_v]);
                cell_labels_2d.push(lab);
            }
        }
    }

    if cells_2d.is_empty() {
        return Err(format!("EME cross-section plane y = {y_cut:.6e} m intersected zero 3D cells"));
    }

    // Intersect 3D facets straddling y = y_cut
    let mut edge_label_map: BTreeMap<[u32; 2], i32> = BTreeMap::new();
    for f in 0..d3.n_facets() {
        let fv = d3.facet(f);
        let mut neg = Vec::with_capacity(3);
        let mut pos = Vec::with_capacity(3);
        for &v in fv {
            if d3.point(v as usize)[1] < y_cut {
                neg.push(v);
            } else {
                pos.push(v);
            }
        }
        let pair = if neg.len() == 1 && pos.len() == 2 {
            Some((
                edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[0], pos[0]),
                edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[0], pos[1]),
            ))
        } else if neg.len() == 2 && pos.len() == 1 {
            Some((
                edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[0], pos[0]),
                edge_cut_vertex(d3, y_cut, &mut edge_to_v2, &mut pts_2d, neg[1], pos[0]),
            ))
        } else {
            None
        };
        if let Some((a, b)) = pair {
            if a != b {
                let key = [a.min(b), a.max(b)];
                edge_label_map.entry(key).or_insert(d3.facet_labels[f]);
            }
        }
    }

    let mut facets_2d = Vec::with_capacity(2 * edge_label_map.len());
    let mut facet_labels_2d = Vec::with_capacity(edge_label_map.len());
    for ([u, w], lab) in edge_label_map {
        facets_2d.extend_from_slice(&[u, w]);
        facet_labels_2d.push(lab);
    }

    let mesh = SimplexMesh::new(2, pts_2d, cells_2d, cell_labels_2d, facets_2d, facet_labels_2d);
    mesh.check()?;
    Ok(mesh)
}

#[inline]
fn orient_ccw_2d(pts: &[f64], a: &mut u32, b: &mut u32, c_v: &mut u32) {
    let (ax, az) = (pts[2 * (*a as usize)], pts[2 * (*a as usize) + 1]);
    let (bx, bz) = (pts[2 * (*b as usize)], pts[2 * (*b as usize) + 1]);
    let (cx, cz) = (pts[2 * (*c_v as usize)], pts[2 * (*c_v as usize) + 1]);
    let det = (bx - ax) * (cz - az) - (bz - az) * (cx - ax);
    if det < 0.0 {
        std::mem::swap(b, c_v);
    }
    // Rotate cyclic order so lowest vertex index comes first (preserves CCW orientation)
    if *b < *a && *b < *c_v {
        let tmp = *a;
        *a = *b;
        *b = *c_v;
        *c_v = tmp;
    } else if *c_v < *a && *c_v < *b {
        let tmp = *c_v;
        *c_v = *b;
        *b = *a;
        *a = tmp;
    }
}

#[inline]
fn dist2_2d(pts: &[f64], u: u32, w: u32) -> f64 {
    let dx = pts[2 * (u as usize)] - pts[2 * (w as usize)];
    let dz = pts[2 * (u as usize) + 1] - pts[2 * (w as usize) + 1];
    dx * dx + dz * dz
}

/// Canonically orders vertices, cells, and facets of a 2D slice so identical cross-sections share exact numbering.
fn canonicalize_2d_mesh(m: &SimplexMesh) -> SimplexMesh {
    let sorted = femcore::meshing::geomdrive::sort_nodes(m);
    let nc = sorted.n_cells();
    let mut cell_rows: Vec<([u32; 3], i32)> = (0..nc)
        .map(|ci| {
            let cv = sorted.cell(ci);
            let (mut a, mut b, mut c_v) = (cv[0], cv[1], cv[2]);
            orient_ccw_2d(&sorted.points, &mut a, &mut b, &mut c_v);
            ([a, b, c_v], sorted.cell_labels[ci])
        })
        .collect();
    cell_rows.sort_by_key(|(v, lab)| (*v, *lab));

    let nf = sorted.n_facets();
    let mut facet_rows: Vec<([u32; 2], i32)> = (0..nf)
        .map(|fi| {
            let fv = sorted.facet(fi);
            ([fv[0].min(fv[1]), fv[0].max(fv[1])], sorted.facet_labels[fi])
        })
        .collect();
    facet_rows.sort_by_key(|(v, lab)| (*v, *lab));

    let mut cells = Vec::with_capacity(3 * nc);
    let mut cell_labels = Vec::with_capacity(nc);
    for (v, lab) in cell_rows {
        cells.extend_from_slice(&v);
        cell_labels.push(lab);
    }
    let mut facets = Vec::with_capacity(2 * nf);
    let mut facet_labels = Vec::with_capacity(nf);
    for (v, lab) in facet_rows {
        facets.extend_from_slice(&v);
        facet_labels.push(lab);
    }
    SimplexMesh::new(2, sorted.points, cells, cell_labels, facets, facet_labels)
}

fn slices_geometrically_compatible(a: &SimplexMesh, b: &SimplexMesh) -> bool {
    if a.n_points() != b.n_points() || a.n_cells() != b.n_cells() || a.cells != b.cells {
        return false;
    }
    let (lo, hi) = a.bbox();
    let diag = (hi[0] - lo[0]).hypot(hi[1] - lo[1]).max(1e-12);
    let tol = 1e-8 * diag;
    a.points.iter().zip(&b.points).all(|(&pa, &pb)| (pa - pb).abs() <= tol)
}

/// Samples domain-filtered temperatures `T(x_q, y_c, z_q)` at Gauss quadrature rules (`quad4` and `quad6`) of `mesh_2d`.
fn sample_slice_temperatures(
    mesh_2d: &SimplexMesh,
    y_c: f64,
    problem: &Problem,
    t_idx: &TransferIndex<'_>,
    t_sf: &SourceField<'_>,
    quad4: &femcore::fem::Quad,
    quad6: &femcore::fem::Quad,
) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>, bool), String> {
    let nc = mesh_2d.n_cells();
    let nq4 = quad4.w.len();
    let nq6 = quad6.w.len();
    let dom_names_cache: BTreeMap<i32, String> =
        mesh_2d.cell_labels.iter().copied().map(|l| (l, canonical_name(&problem.domain_names, l))).collect();

    let eval_rule = |qrule: &femcore::fem::Quad| -> Result<Vec<Vec<f64>>, String> {
        let nq = qrule.w.len();
        let mut coords = Vec::with_capacity(nc * nq * 3);
        let mut filters: Vec<Option<&str>> = Vec::with_capacity(nc * nq);
        for ci in 0..nc {
            let vs = mesh_2d.cell(ci);
            let (p0, p1, p2) =
                (mesh_2d.point(vs[0] as usize), mesh_2d.point(vs[1] as usize), mesh_2d.point(vs[2] as usize));
            let dname = dom_names_cache.get(&mesh_2d.cell_labels[ci]).map(String::as_str);
            for l in &qrule.lam {
                let x = l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0];
                let z = l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1];
                coords.extend_from_slice(&[x, y_c, z]);
                filters.push(dname);
            }
        }
        let sampled = t_idx.eval_at_points(t_sf, &coords, &filters).map_err(|e| e.to_string())?;
        let mut out = Vec::with_capacity(nc);
        for ci in 0..nc {
            let dom = mesh_2d.cell_labels[ci];
            let fallback = problem
                .fallback_temperature
                .or_else(|| problem.thermo_optic_domains.get(&dom).and_then(|td| td.spec.reference_temperature))
                .unwrap_or(293.15);
            let mut row = Vec::with_capacity(nq);
            for q in 0..nq {
                let idx = ci * nq + q;
                let tk = if sampled.found[idx] { sampled.values[idx] } else { fallback };
                row.push(tk);
            }
            out.push(row);
        }
        Ok(out)
    };

    let q4 = eval_rule(quad4)?;
    let q6 = eval_rule(quad6)?;
    let mut avgs = Vec::with_capacity(nc);
    let mut subcell_grad = false;
    for ci in 0..nc {
        let mut mean: f64 = 0.0;
        for q in 0..nq4 {
            mean += quad4.w[q] * q4[ci][q];
        }
        // note: on 2D reference triangle, sum(quad4.w) == 1.0 in femcore (area is multiplied separately)
        let wsum: f64 = quad4.w.iter().sum::<f64>().max(1e-30);
        mean /= wsum;
        for q in 0..nq6 {
            if (q6[ci][q] - mean).abs() > 1e-6 * mean.abs().max(1.0) {
                subcell_grad = true;
            }
        }
        avgs.push(mean);
    }
    Ok((q4, q6, avgs, subcell_grad))
}

enum VecSliceSource<'a, const K: usize> {
    Nodal { sf: SourceField<'a>, src_dim: usize, ncomp: usize },
    PerCell3d(&'a [[f64; K]]),
}

/// Samples a `K`-component vector/Voigt field (`K = 3` for wafer `E_DC`, `K = 6` for solver-Voigt `S`)
/// at `quad4` and `quad6` of `mesh_2d` with strict canonical `domain_names` filtering.
fn sample_slice_vec_field<const K: usize>(
    mesh_2d: &SimplexMesh,
    y_c: f64,
    problem: &Problem,
    idx: &TransferIndex<'_>,
    src: &VecSliceSource<'_, K>,
    quad4: &femcore::fem::Quad,
    quad6: &femcore::fem::Quad,
) -> Result<(Vec<Vec<[f64; K]>>, Vec<Vec<[f64; K]>>, Vec<[f64; K]>, bool), String> {
    let nc = mesh_2d.n_cells();
    let nq4 = quad4.w.len();
    let nq6 = quad6.w.len();
    let dom_names_cache: BTreeMap<i32, String> =
        mesh_2d.cell_labels.iter().copied().map(|l| (l, canonical_name(&problem.domain_names, l))).collect();

    let eval_rule = |qrule: &femcore::fem::Quad| -> Result<Vec<Vec<[f64; K]>>, String> {
        let nq = qrule.w.len();
        match src {
            VecSliceSource::Nodal { sf, src_dim, ncomp } => {
                let step_dim = if *src_dim == 2 { 2 } else { 3 };
                let mut coords = Vec::with_capacity(nc * nq * step_dim);
                let mut filters: Vec<Option<&str>> = Vec::with_capacity(nc * nq);
                for ci in 0..nc {
                    let vs = mesh_2d.cell(ci);
                    let (p0, p1, p2) =
                        (mesh_2d.point(vs[0] as usize), mesh_2d.point(vs[1] as usize), mesh_2d.point(vs[2] as usize));
                    let dname = dom_names_cache.get(&mesh_2d.cell_labels[ci]).map(String::as_str);
                    for l in &qrule.lam {
                        let x = l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0];
                        let z = l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1];
                        if *src_dim == 2 {
                            coords.extend_from_slice(&[x, z]);
                        } else {
                            coords.extend_from_slice(&[x, y_c, z]);
                        }
                        filters.push(dname);
                    }
                }
                let sampled = idx.eval_at_points(sf, &coords, &filters).map_err(|e| e.to_string())?;
                let mut out = Vec::with_capacity(nc);
                for ci in 0..nc {
                    let mut row = Vec::with_capacity(nq);
                    for q in 0..nq {
                        let p_idx = ci * nq + q;
                        let mut val = [0.0; K];
                        if sampled.found[p_idx] {
                            let base = p_idx * (*ncomp);
                            if K == 3 && *src_dim == 2 && *ncomp == 2 {
                                val[0] = sampled.values[base];
                                val[2] = sampled.values[base + 1];
                            } else if *ncomp >= K {
                                val[..K].copy_from_slice(&sampled.values[base..(base + K)]);
                            }
                        }
                        row.push(val);
                    }
                    out.push(row);
                }
                Ok(out)
            }
            VecSliceSource::PerCell3d(cells) => {
                let mut out = Vec::with_capacity(nc);
                for ci in 0..nc {
                    let vs = mesh_2d.cell(ci);
                    let (p0, p1, p2) =
                        (mesh_2d.point(vs[0] as usize), mesh_2d.point(vs[1] as usize), mesh_2d.point(vs[2] as usize));
                    let dname = dom_names_cache.get(&mesh_2d.cell_labels[ci]).map(String::as_str);
                    let mut row = Vec::with_capacity(nq);
                    for l in &qrule.lam {
                        let x = l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0];
                        let z = l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1];
                        let val = match idx.locate(&[x, y_c, z], dname) {
                            Some(loc) if loc.cell < cells.len() => cells[loc.cell],
                            _ => [0.0; K],
                        };
                        row.push(val);
                    }
                    out.push(row);
                }
                Ok(out)
            }
        }
    };

    let q4 = eval_rule(quad4)?;
    let q6 = eval_rule(quad6)?;
    let wsum: f64 = quad4.w.iter().sum::<f64>().max(1e-30);
    let mut avgs = Vec::with_capacity(nc);
    let mut subcell_grad = false;
    for ci in 0..nc {
        let mut mean = [0.0; K];
        for q in 0..nq4 {
            for a in 0..K {
                mean[a] += quad4.w[q] * q4[ci][q][a];
            }
        }
        let mut scale_ref = 0.0_f64;
        for a in 0..K {
            mean[a] /= wsum;
            scale_ref = scale_ref.max(mean[a].abs());
        }
        for q in 0..nq6 {
            for a in 0..K {
                scale_ref = scale_ref.max(q6[ci][q][a].abs());
            }
        }
        if scale_ref > 0.0 {
            for q in 0..nq6 {
                for a in 0..K {
                    if (q6[ci][q][a] - mean[a]).abs() > 1e-6 * scale_ref {
                        subcell_grad = true;
                    }
                }
            }
        }
        avgs.push(mean);
    }
    Ok((q4, q6, avgs, subcell_grad))
}

#[inline]
fn apply_eo_pe_to_c3(
    mut ep: C3,
    r_g: Option<&[[f64; 3]; 6]>,
    p_g: Option<&femcore::material::M6>,
    e_dc: Option<&[f64; 3]>,
    s_dc: Option<&[f64; 6]>,
) -> C3 {
    if (r_g.is_none() || e_dc.is_none()) && (p_g.is_none() || s_dc.is_none()) {
        return ep;
    }
    let re_cur = [
        [ep[0][0].re, ep[0][1].re, ep[0][2].re],
        [ep[1][0].re, ep[1][1].re, ep[1][2].re],
        [ep[2][0].re, ep[2][1].re, ep[2][2].re],
    ];
    let d_eps = femcore::material::coupled::delta_eps_eo_pe(&re_cur, r_g, p_g, e_dc, s_dc);
    for a in 0..3 {
        for b in 0..3 {
            ep[a][b] += c(d_eps[a][b], 0.0);
        }
    }
    ep
}

/// Samples electrostatic potential `V(x, y_c, z)` at the 2D P2 Lagrange DOF coordinates of `mesh_2d`.
fn sample_slice_potential(
    mesh_2d: &SimplexMesh,
    y_c: f64,
    p_idx: &TransferIndex<'_>,
    p_sf: &SourceField<'_>,
) -> Result<Vec<f64>, String> {
    let coords_2d = mesh_2d.dof_coordinates(2);
    let nd = coords_2d.len() / 2;
    let mut coords_3d = Vec::with_capacity(nd * 3);
    for i in 0..nd {
        coords_3d.extend_from_slice(&[coords_2d[2 * i], y_c, coords_2d[2 * i + 1]]);
    }
    let sampled = p_idx.eval_at_points(p_sf, &coords_3d, &[]).map_err(|e| e.to_string())?;
    Ok(sampled.values)
}

/// Computes a deterministic BLAKE3 digest of a 2D slice's canonical geometry, boundary conditions,
/// cylindrical bend radius, local constitutive tensors (including thermo-optic `T`, electro-optic `E_DC`,
/// and photoelastic `S`), and gated sheets.
#[allow(clippy::too_many_arguments)]
fn compute_slice_blake3(
    mesh_2d: &SimplexMesh,
    problem: &Problem,
    freq: f64,
    bend_radius: Option<f64>,
    eo_pe_domains: &BTreeMap<i32, (Option<[[f64; 3]; 6]>, Option<femcore::material::M6>)>,
    cell_quad4_temps: &[Vec<f64>],
    cell_quad4_e_dc: &[Vec<[f64; 3]>],
    cell_quad4_strain: &[Vec<[f64; 6]>],
    slice_potential: Option<&[f64]>,
) -> Result<String, String> {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"nyx-eme-slice-v1");
    hasher.update(&(problem.order as u64).to_le_bytes());
    hasher.update(&quantize_f64(freq, 1e-9).to_le_bytes());
    hasher.update(problem.symmetry.as_str().as_bytes());

    match bend_radius {
        Some(r) => {
            hasher.update(&[1u8]);
            hasher.update(&quantize_f64(r, 1e-12).to_le_bytes());
        }
        None => {
            hasher.update(&[0u8]);
        }
    }

    let (lo, hi) = mesh_2d.bbox();
    let diag = (hi[0] - lo[0]).hypot(hi[1] - lo[1]).max(1e-12);
    let coord_unit = 1e-11 * diag;
    hasher.update(&(mesh_2d.n_points() as u64).to_le_bytes());
    for &c_val in &mesh_2d.points {
        hasher.update(&quantize_f64(c_val, coord_unit).to_le_bytes());
    }
    hasher.update(&(mesh_2d.cells.len() as u64).to_le_bytes());
    for &vi in &mesh_2d.cells {
        hasher.update(&vi.to_le_bytes());
    }
    for &cl in &mesh_2d.cell_labels {
        hasher.update(&cl.to_le_bytes());
    }
    hasher.update(&(mesh_2d.facets.len() as u64).to_le_bytes());
    for &fi in &mesh_2d.facets {
        hasher.update(&fi.to_le_bytes());
    }
    for &fl in &mesh_2d.facet_labels {
        hasher.update(&fl.to_le_bytes());
    }

    // Hash evaluated constitutive tensors per cell (incorporating local thermo-optic, EO, and PE perturbations!)
    let eps_base = problem.materials.cell_full_complex_permittivity(mesh_2d, Some(1.0))?;
    let mui_base = problem.materials.cell_full_complex_permeability_inv(mesh_2d, Some(1.0))?;
    for ci in 0..mesh_2d.n_cells() {
        let dom = mesh_2d.cell_labels[ci];
        let td_opt = problem.thermo_optic_domains.get(&dom);
        let eo_pe_opt = eo_pe_domains.get(&dom);
        let qtemps = cell_quad4_temps.get(ci);
        let qe = cell_quad4_e_dc.get(ci);
        let qs = cell_quad4_strain.get(ci);
        let has_to = td_opt.is_some() && qtemps.is_some();
        let has_eope = eo_pe_opt.is_some() && (qe.is_some() || qs.is_some());
        if has_to || has_eope {
            let nq = qtemps.map_or(0, Vec::len).max(qe.map_or(0, Vec::len)).max(qs.map_or(0, Vec::len));
            for q in 0..nq {
                let mut ep = eps_base[ci];
                if let (Some(td), Some(qt)) = (td_opt, qtemps) {
                    if let Some(&tk) = qt.get(q) {
                        if tk > 0.0 && tk.is_finite() {
                            if let Some((re, im)) = td.eval_at(tk)? {
                                ep = apply_thermo_optic_to_c3(ep, &re, &im);
                            }
                        }
                    }
                }
                if let Some((r_g, p_g)) = eo_pe_opt {
                    let e_q = qe.and_then(|v| v.get(q));
                    let s_q = qs.and_then(|v| v.get(q));
                    ep = apply_eo_pe_to_c3(ep, r_g.as_ref(), p_g.as_ref(), e_q, s_q);
                }
                hash_c3(&mut hasher, &ep);
            }
        } else {
            hash_c3(&mut hasher, &eps_base[ci]);
        }
        hash_c3(&mut hasher, &mui_base[ci]);
        if let Some(pml) = problem.materials.pml.get(&dom) {
            hasher.update(format!("{pml:?}").as_bytes());
        }
    }

    for &pl in &problem.pec {
        hasher.update(&pl.to_le_bytes());
    }
    for &pml_l in &problem.pmc {
        hasher.update(&pml_l.to_le_bytes());
    }
    for ibc in &problem.surface_impedance {
        hasher.update(format!("{ibc:?}").as_bytes());
    }
    for sh in &problem.sheet_2d {
        hasher.update(format!("{sh:?}").as_bytes());
    }
    if let Some(pot) = slice_potential {
        hasher.update(&(pot.len() as u64).to_le_bytes());
        for &pv in pot {
            hasher.update(&quantize_f64(pv, 1e-10).to_le_bytes());
        }
    }

    Ok(hasher.finalize().to_hex().to_string())
}

#[inline]
fn quantize_f64(v: f64, scale: f64) -> i64 {
    let unit = scale.abs().max(1e-30);
    let q = (v / unit).round();
    if q >= (i64::MAX as f64) {
        i64::MAX
    } else if q <= (i64::MIN as f64) {
        i64::MIN
    } else {
        q as i64
    }
}

#[inline]
fn hash_c3(hasher: &mut blake3::Hasher, m: &C3) {
    for r in 0..3 {
        for col in 0..3 {
            hasher.update(&quantize_f64(m[r][col].re, 1e-12).to_le_bytes());
            hasher.update(&quantize_f64(m[r][col].im, 1e-12).to_le_bytes());
        }
    }
}

/// Extracts 2D `SimplexMesh` patches for any explicit `waveport` BCs lying on `y = y0` or `y = y1`.
fn extract_eme_waveport_patches(
    d3: &SimplexMesh,
    waveports: &[WavePortBc],
    pmc_labels: &[i32],
    y0: f64,
    y1: f64,
) -> Result<Vec<EmeWavePortPatch>, String> {
    if waveports.is_empty() {
        return Ok(Vec::new());
    }
    let span_y = (y1 - y0).abs().max(1e-12);
    let tol_y = 1e-6 * span_y;
    let face_to_cell = build_exterior_face_to_cell(d3);
    let exterior = d3.exterior_facets();
    let mut pmc_edges_3d = BTreeSet::new();
    for f in 0..d3.n_facets() {
        if pmc_labels.contains(&d3.facet_labels[f]) {
            let fv = d3.facet(f);
            for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[2], fv[0])] {
                pmc_edges_3d.insert((a.min(b), a.max(b)));
            }
        }
    }
    let mut out = Vec::with_capacity(waveports.len());

    for wp in waveports {
        let mut tris = Vec::new();
        let mut on_y0 = true;
        let mut on_y1 = true;
        for f in 0..d3.n_facets() {
            if !wp.labels.contains(&d3.facet_labels[f]) || !exterior[f] {
                continue;
            }
            let fv = d3.facet(f);
            for &v in fv {
                let y = d3.point(v as usize)[1];
                if (y - y0).abs() > tol_y {
                    on_y0 = false;
                }
                if (y - y1).abs() > tol_y {
                    on_y1 = false;
                }
            }
            tris.push(f);
        }
        if tris.is_empty() {
            return Err(format!("EME waveport '{}' selected zero exterior boundary triangles", wp.name));
        }
        if !on_y0 && !on_y1 {
            return Err(format!(
                "EME waveport '{}' must lie entirely on the input plane y = {y0:.6e} m or output plane y = {y1:.6e} m",
                wp.name
            ));
        }
        let is_output = on_y1 && !on_y0;

        let mut vmap: HashMap<u32, u32> = HashMap::new();
        let mut v2_to_v3: Vec<u32> = Vec::new();
        let mut pts_2d = Vec::new();
        let mut cells_2d = Vec::with_capacity(3 * tris.len());
        let mut cell_labels_2d = Vec::with_capacity(tris.len());
        let mut edge_counts: BTreeMap<[u32; 2], usize> = BTreeMap::new();

        for &f in &tris {
            let fv = d3.facet(f);
            let mut lv = [0u32; 3];
            for k in 0..3 {
                let gv = fv[k];
                lv[k] = *vmap.entry(gv).or_insert_with(|| {
                    let p = d3.point(gv as usize);
                    let id = (pts_2d.len() / 2) as u32;
                    pts_2d.push(p[0]);
                    pts_2d.push(p[2]);
                    v2_to_v3.push(gv);
                    id
                });
            }
            let (mut a, mut b, mut c_v) = (lv[0], lv[1], lv[2]);
            orient_ccw_2d(&pts_2d, &mut a, &mut b, &mut c_v);
            cells_2d.extend_from_slice(&[a, b, c_v]);
            let key3 = sorted_tri_key(fv[0], fv[1], fv[2]);
            let ci = face_to_cell.get(&key3).copied().unwrap_or(0);
            cell_labels_2d.push(d3.cell_labels[ci]);

            for (u, w) in [(a, b), (b, c_v), (c_v, a)] {
                *edge_counts.entry([u.min(w), u.max(w)]).or_insert(0) += 1;
            }
        }

        // Exterior rim edges of the waveport patch get PMC (label 2) if on a 3D PMC wall, else PEC (label 1)
        let mut facets_2d = Vec::new();
        let mut facet_labels_2d = Vec::new();
        for ([u, w], cnt) in edge_counts {
            if cnt == 1 {
                let (ga, gb) = (v2_to_v3[u as usize], v2_to_v3[w as usize]);
                let e3 = (ga.min(gb), ga.max(gb));
                facets_2d.extend_from_slice(&[u, w]);
                facet_labels_2d.push(if pmc_edges_3d.contains(&e3) { 2 } else { 1 });
            }
        }
        let mesh_2d =
            canonicalize_2d_mesh(&SimplexMesh::new(2, pts_2d, cells_2d, cell_labels_2d, facets_2d, facet_labels_2d));
        mesh_2d.check()?;
        out.push(EmeWavePortPatch { bc: wp.clone(), is_output, mesh_2d });
    }
    Ok(out)
}

fn sorted_tri_key(a: u32, b: u32, c_v: u32) -> [u32; 3] {
    let mut t = [a, b, c_v];
    t.sort_unstable();
    t
}

fn build_exterior_face_to_cell(d3: &SimplexMesh) -> HashMap<[u32; 3], usize> {
    let mut map = HashMap::new();
    let faces = [[1usize, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    for ci in 0..d3.n_cells() {
        let cv = d3.cell(ci);
        for f in &faces {
            map.entry(sorted_tri_key(cv[f[0]], cv[f[1]], cv[f[2]])).or_insert(ci);
        }
    }
    map
}

// ------------------------------------------------------------------------------------------------------------- Stage 2: Slice Modes, Overlaps & Cascade

/// Stage 2: solves the unique 2D cross-section eigenproblems, computes exact cross-power interface overlaps and
/// `S_if` matrices, and evaluates the nominal multi-port scattering cascade via Redheffer star products.
pub fn solve_eme_cascade(asm: &EmeAssembly) -> Result<EmeCascade, String> {
    let n_unique = asm.unique_slices.len();
    let mut slice_modes = Vec::with_capacity(n_unique);
    let mut max_residual = 0.0_f64;
    let mut total_restarts = 0usize;
    let mut all_converged = true;

    for us in &asm.unique_slices {
        femcore::progress::check()?;
        let solved = solve_unique_slice_modes(us, &asm.problem, asm.freq_hz, asm.n_target, asm.evanescent)?;
        for &r in &solved.residuals {
            max_residual = max_residual.max(r);
        }
        total_restarts += solved.restarts;
        all_converged &= solved.converged;
        slice_modes.push(solved);
    }

    // Collect every ordered (left_slice, right_slice) interface needed by the cascade tree (including Periodic wraps)
    let mut needed_pairs: Vec<(usize, usize)> = Vec::new();
    collect_tree_interface_pairs(&asm.tree, &mut needed_pairs);
    needed_pairs.sort_unstable();
    needed_pairs.dedup();

    let mut interfaces = Vec::with_capacity(needed_pairs.len());
    for (u_l, u_r) in needed_pairs {
        femcore::progress::check()?;
        let ov = compute_interface_overlap(u_l, &slice_modes[u_l], u_r, &slice_modes[u_r])?;
        interfaces.push(ov);
    }

    // Optional explicit waveports on y = y0 and y = y1
    let (
        input_port_modes,
        output_port_modes,
        input_port_overlap,
        output_port_overlap,
        port_channels,
        ports,
        port_z0_ohm,
    ) = if !asm.waveport_patches.is_empty() {
        let (first_u, last_u) = tree_end_slices(&asm.tree);
        build_explicit_waveport_couplings(asm, &slice_modes[first_u], first_u, &slice_modes[last_u], last_u)?
    } else {
        let (first_u, last_u) = tree_end_slices(&asm.tree);
        let m_in = slice_modes[first_u].neff_fwd.len();
        let m_out = slice_modes[last_u].neff_fwd.len();
        let p_in = if asm.port_modes == 0 { m_in } else { asm.port_modes.min(m_in) };
        let p_out = if asm.port_modes == 0 { m_out } else { asm.port_modes.min(m_out) };
        let mut ports = Vec::with_capacity(p_in + p_out);
        let mut channels = Vec::with_capacity(p_in + p_out);
        for m in 0..p_in {
            ports.push(format!("p1_m{}", m + 1));
            channels.push((false, m));
        }
        for m in 0..p_out {
            ports.push(format!("p2_m{}", m + 1));
            channels.push((true, m));
        }
        let z0 = vec![asm.port_z0_ohm; ports.len()];
        (None, None, None, None, channels, ports, z0)
    };

    let summaries: Vec<EmeUniqueSliceSummary> = asm
        .unique_slices
        .iter()
        .zip(&slice_modes)
        .map(|(us, sm)| EmeUniqueSliceSummary {
            blake3_hash: us.blake3_hash.clone(),
            sample_y_m: us.sample_y_m,
            n_modes: sm.neff_fwd.len(),
            bend_radius_m: us.bend_radius_m,
            neff_fwd: sm.neff_fwd.iter().map(|z| [z.re, z.im]).collect(),
            neff_bwd: sm.neff_bwd.iter().map(|z| [z.re, z.im]).collect(),
            beta_fwd: sm.beta_fwd.iter().map(|z| [z.re, z.im]).collect(),
            beta_bwd: sm.beta_bwd.iter().map(|z| [z.re, z.im]).collect(),
            d_self: sm.d_self.iter().map(|z| [z.re, z.im]).collect(),
            d_self_tilde: sm.d_self_tilde.iter().map(|z| [z.re, z.im]).collect(),
        })
        .collect();

    let data = EmeData {
        freq_hz: asm.freq_hz,
        y_range_m: asm.y_range_m,
        total_segments: asm.mesh_segments.len(),
        total_periodic_subgrids: asm.total_periodic_subgrids,
        unique_slices: summaries,
        tree: asm.tree.clone(),
        interfaces,
        input_port_overlap,
        output_port_overlap,
        port_channels,
        ports,
        port_z0_ohm,
    };

    let s_terminal = evaluate_eme_terminal_s(&data, None, None, None)?;

    Ok(EmeCascade {
        slice_modes,
        input_port_modes,
        output_port_modes,
        data,
        s_terminal,
        max_residual,
        total_restarts,
        all_converged,
    })
}

fn solve_unique_slice_modes(
    us: &EmeUniqueSlice,
    problem: &Problem,
    freq: f64,
    n_target: f64,
    evanescent: bool,
) -> Result<SolvedSliceModes, String> {
    let mesh = &us.mesh_2d;
    let order = problem.order;
    let nc = mesh.n_cells();
    let k0 = 2.0 * PI * freq / C0;

    let pml_table = EmPmlTable::from_simplex(mesh, &problem.materials, freq, true, false);
    let has_pml = pml_table.has_any();
    let pml_cells: Vec<bool> =
        (0..nc).map(|ci| problem.materials.pml.get(&mesh.cell_labels[ci]).is_some_and(Pml::is_active)).collect();

    let eps_base_w = problem.materials.cell_full_complex_permittivity(mesh, Some(1.0))?;
    let mui_base_w = problem.materials.cell_full_complex_permeability_inv(mesh, Some(1.0))?;

    let eo_pe_domains = problem.domain_eo_pe_global();

    // Check whether the fast real-symmetric waveguide eigensolver applies to this slice
    let can_try_real = !has_pml
        && us.bend_radius_m.is_none()
        && !us.has_subcell_temp_grad
        && !us.has_subcell_eo_pe_grad
        && problem.sheet_2d.is_empty()
        && !problem.has_loss()
        && !problem.materials.is_non_symmetric()
        && !is_c3_slice_non_symmetric(&eps_base_w, &mui_base_w)
        && problem.materials.aniso.values().all(|a| {
            a.permittivity.as_ref().is_none_or(|t| {
                let loc = local_tensor(&femcore::material::tensor_global(t, a.euler_deg));
                loc[0][2].abs() <= 1e-14 && loc[1][2].abs() <= 1e-14
            })
        });

    let mut mui_phys_wafer = mui_base_w.clone();
    if let Some(r_bend) = us.bend_radius_m {
        for ci in 0..nc {
            let cv = mesh.cell(ci);
            let xc =
                (mesh.point(cv[0] as usize)[0] + mesh.point(cv[1] as usize)[0] + mesh.point(cv[2] as usize)[0]) / 3.0;
            let hy = (1.0 + xc / r_bend).max(1e-3);
            mui_phys_wafer[ci] = apply_bend_metric_inv_wafer(&mui_phys_wafer[ci], hy);
        }
    }

    if can_try_real {
        // Build cellwise real EpsCell (incorporating cell-average thermo-optic, EO, and PE perturbations when present)
        let mut eps_cells = Vec::with_capacity(nc);
        let mut real_ok = true;
        for ci in 0..nc {
            let mut ep_w = eps_base_w[ci];
            let dom = mesh.cell_labels[ci];
            if let (Some(td), Some(&tk)) = (problem.thermo_optic_domains.get(&dom), us.cell_avg_temps.get(ci)) {
                if tk > 0.0 && tk.is_finite() {
                    if let Some((re, im)) = td.eval_at(tk)? {
                        ep_w = apply_thermo_optic_to_c3(ep_w, &re, &im);
                    }
                }
            }
            if let Some((r_g, p_g)) = eo_pe_domains.get(&dom) {
                ep_w = apply_eo_pe_to_c3(
                    ep_w,
                    r_g.as_ref(),
                    p_g.as_ref(),
                    us.cell_avg_e_dc.get(ci),
                    us.cell_avg_strain.get(ci),
                );
            }
            let re_w = [
                [ep_w[0][0].re, ep_w[0][1].re, ep_w[0][2].re],
                [ep_w[1][0].re, ep_w[1][1].re, ep_w[1][2].re],
                [ep_w[2][0].re, ep_w[2][1].re, ep_w[2][2].re],
            ];
            let has_im = ep_w.iter().flatten().any(|z| z.im.abs() > 1e-14);
            let loc = local_tensor(&re_w);
            if has_im || loc[0][2].abs() > 1e-14 || loc[1][2].abs() > 1e-14 {
                real_ok = false;
                break;
            }
            eps_cells.push(EpsCell { xx: loc[0][0], xy: loc[0][1], yy: loc[1][1], zz: loc[2][2] });
        }
        if real_ok {
            let mur = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
            let r_asm = assemble_modes_with_bcs_and_symmetry(
                mesh,
                freq,
                &eps_cells,
                &mur,
                order,
                &problem.pec,
                &problem.surface_impedance,
                &problem.pmc,
                problem.symmetry,
            )?;
            let k_req = us.max_modes.min(r_asm.a.nrows().saturating_sub(1)).max(1);
            if let Ok(r_modes) = solve_mode_pairs(&r_asm, n_target, evanescent, k_req) {
                if r_modes.converged && !r_modes.complex_wanted && r_modes.neff.len() >= k_req {
                    let m_got = r_modes.neff.len().min(us.max_modes);
                    let neff_c: Vec<C> = (0..m_got).map(|i| c(r_modes.neff[i].0, r_modes.neff[i].1)).collect();
                    let et_c: Vec<Vec<C>> =
                        (0..m_got).map(|i| r_modes.et[i].iter().map(|&v| c(v, 0.0)).collect()).collect();
                    let phi_c: Vec<Vec<C>> =
                        (0..m_got).map(|i| r_modes.phi[i].iter().map(|&v| c(v, 0.0)).collect()).collect();
                    let mc = ModesC {
                        neff: neff_c,
                        et: et_c.clone(),
                        phi: phi_c.clone(),
                        et_left: et_c,
                        phi_left: phi_c,
                        biorth_norms: vec![c(1.0, 0.0); m_got],
                        is_non_symmetric: false,
                        residuals: r_modes.residuals[..m_got].to_vec(),
                        restarts: r_modes.restarts,
                        matvecs: r_modes.matvecs,
                        converged: r_modes.converged,
                    };
                    return finalize_slice_modes(
                        us,
                        problem,
                        freq,
                        k0,
                        false,
                        &pml_cells,
                        &eps_base_w,
                        &mui_phys_wafer,
                        mc,
                        None,
                    );
                }
            }
        }
    }

    // Full quadrature-point complex waveguide mode assembly (PML, bend_radius, subcell T(x,z), gyrotropy, sheet_2d, loss)
    let c_asm = assemble_slice_complex_quad(us, problem, freq, &pml_table, &eps_base_w, &mui_base_w)?;
    let n_free = c_asm.a.nrows();
    if n_free < 2 {
        return Err("EME 2D cross-section has fewer than 2 free electromagnetic DOFs".into());
    }
    let k_ask = (us.max_modes + 4).min(n_free - 1).max(1);
    let target = if evanescent { c(0.0, n_target) } else { c(n_target, 0.0) };
    let mut mc_fwd = solve_mode_pairs_complex(&c_asm, target, k_ask)?;
    if mc_fwd.neff.len() < us.max_modes && k_ask < n_free - 1 {
        let k_retry = (2 * us.max_modes + 6).min(n_free - 1);
        if let Ok(mc_more) = solve_mode_pairs_complex(&c_asm, target, k_retry) {
            if mc_more.neff.len() > mc_fwd.neff.len() {
                mc_fwd = mc_more;
            }
        }
    }
    if mc_fwd.neff.is_empty() {
        return Err("EME 2D complex eigensolver returned zero converged modes on cross-section slice".into());
    }

    // For non-reciprocal Voigt / magneto-optic waveguides with longitudinal-transverse A1 coupling,
    // also solve backward (-n_target) propagation constants so beta_bwd != beta_fwd is captured!
    let mc_bwd = if c_asm.is_non_symmetric && c_asm.a1.is_some() && !evanescent {
        solve_mode_pairs_complex(&c_asm, -target, k_ask).ok()
    } else {
        None
    };

    finalize_slice_modes(us, problem, freq, k0, has_pml, &pml_cells, &eps_base_w, &mui_phys_wafer, mc_fwd, mc_bwd)
}

/// Quadrature-point complex 2D waveguide assembly incorporating `ThermoOpticDomain` `T(x_q, y_c, z_q)`,
/// Ward–Pendry cylindrical bend metric `h_y(x_q) = 1 + x_q / R_bend`, `EmPmlTable`, `surface_impedance`,
/// `sheet_2d` (+ `slice_potential`), `pec`, `pmc`, and midplane `symmetry`.
fn assemble_slice_complex_quad(
    us: &EmeUniqueSlice,
    problem: &Problem,
    freq: f64,
    pml: &EmPmlTable,
    eps_base_w: &[C3],
    mui_base_w: &[C3],
) -> Result<ModeAssemblyC, String> {
    let mesh = &us.mesh_2d;
    let order = problem.order;
    let (nc, nv, ne) = (mesh.n_cells(), mesh.n_points(), mesh.n_edges());
    let k0 = 2.0 * PI * freq / C0;
    let k2 = c(k0 * k0, 0.0);
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let net = cd.net;
    let nphi = mesh.n_dofs(order);
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let quad4 = quad_rule(2, 4);
    let quad6 = quad_rule(2, 6);
    let tab4: Vec<(Vec<f64>, Vec<f64>)> = quad4.lam.iter().map(|l| shape(l, order)).collect();
    let tab6: Vec<(Vec<f64>, Vec<f64>)> = quad6.lam.iter().map(|l| shape(l, order)).collect();

    let perm_loc = |t: &C3| -> C3 {
        let axes = [0usize, 2, 1];
        let mut r = [[c(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = t[axes[i]][axes[j]];
            }
        }
        r
    };

    let is_non_symmetric = problem.materials.is_non_symmetric()
        || is_c3_slice_non_symmetric(eps_base_w, mui_base_w)
        || problem.sheet_2d.iter().any(Sheet2dBc::is_non_symmetric);

    let n_tot = net + nphi;
    let mut a_coo = CCoo::new(n_tot, n_tot);
    let mut a1_coo = CCoo::new(n_tot, n_tot);
    let mut b_coo = CCoo::new(n_tot, n_tot);
    let mut has_coupling = false;

    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut gn = vec![[0.0; 2]; nl];

    let eo_pe_domains = problem.domain_eo_pe_global();
    let use_q6_all = us.bend_radius_m.is_some() || us.has_subcell_temp_grad || us.has_subcell_eo_pe_grad;

    for ci in 0..nc {
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        let is_pml = pml.is_pml(ci);
        let dom = mesh.cell_labels[ci];
        let vs = mesh.cell(ci);
        let pts = [mesh.point(vs[0] as usize), mesh.point(vs[1] as usize), mesh.point(vs[2] as usize)];
        let use_q6 = is_pml || use_q6_all;
        let (qrule, tab) = if use_q6 { (&quad6, &tab6) } else { (&quad4, &tab4) };
        let qtemps = if use_q6 { us.cell_quad6_temps.get(ci) } else { us.cell_quad4_temps.get(ci) };
        let qe_dc = if use_q6 { us.cell_quad6_e_dc.get(ci) } else { us.cell_quad4_e_dc.get(ci) };
        let qstrain = if use_q6 { us.cell_quad6_strain.get(ci) } else { us.cell_quad4_strain.get(ci) };
        let td_opt = problem.thermo_optic_domains.get(&dom);
        let eo_pe_opt = eo_pe_domains.get(&dom);

        let mut amat = vec![c(0.0, 0.0); nloc_e * nloc_e];
        let mut tm = vec![c(0.0, 0.0); nloc_e * nloc_e];
        let mut gm = vec![c(0.0, 0.0); nloc_e * nl];
        let mut gt = vec![c(0.0, 0.0); nl * nloc_e];
        let mut pp = vec![c(0.0, 0.0); nl * nl];
        let mut a1_ee = vec![c(0.0, 0.0); nloc_e * nloc_e];
        let mut mg = vec![c(0.0, 0.0); nloc_e * nl];
        let mut mgt = vec![c(0.0, 0.0); nl * nloc_e];
        let mut cell_cpl = false;
        let mut cell_mu_cpl = false;

        for (q, l) in qrule.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
            let (nv_, dn) = &tab[q];
            for a in 0..nl {
                gn[a][0] = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
                gn[a][1] = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
            }
            let xq = [
                l[0] * pts[0][0] + l[1] * pts[1][0] + l[2] * pts[2][0],
                l[0] * pts[0][1] + l[1] * pts[1][1] + l[2] * pts[2][1],
            ];

            let mut ep_w = eps_base_w[ci];
            let mut mui_w = mui_base_w[ci];
            if let (Some(td), Some(qt)) = (td_opt, qtemps) {
                let tk = qt[q];
                if tk > 0.0 && tk.is_finite() {
                    if let Some((re, im)) = td.eval_at(tk)? {
                        ep_w = apply_thermo_optic_to_c3(ep_w, &re, &im);
                    }
                }
            }
            if let Some((r_g, p_g)) = eo_pe_opt {
                ep_w = apply_eo_pe_to_c3(
                    ep_w,
                    r_g.as_ref(),
                    p_g.as_ref(),
                    qe_dc.and_then(|v| v.get(q)),
                    qstrain.and_then(|v| v.get(q)),
                );
            }
            if let Some(r_bend) = us.bend_radius_m {
                let hy = (1.0 + xq[0] / r_bend).max(1e-3);
                ep_w = apply_bend_metric_wafer(&ep_w, hy);
                mui_w = apply_bend_metric_inv_wafer(&mui_w, hy);
            }

            let ep_loc = perm_loc(&ep_w);
            let mui_loc = perm_loc(&mui_w);
            let (ep_st, mui_st, _) = pml.eval_2d(ci, ep_loc, mui_loc, xq, false);

            add_waveguide_complex_point_full(
                nloc_e,
                nl,
                &w,
                &curl,
                nv_,
                &gn,
                qrule.w[q] * area,
                k2,
                &ep_st,
                &mui_st,
                &mut amat,
                &mut tm,
                &mut gm,
                &mut gt,
                &mut pp,
                &mut a1_ee,
                &mut mg,
                &mut mgt,
                &mut cell_cpl,
                &mut cell_mu_cpl,
            );
        }

        mesh.cell_dofs(ci, order, &mut pd);
        let rows_e = &cd.gdofs[ci];
        let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();
        a_coo.push_block(rows_e, rows_e, &amat);
        b_coo.push_block(rows_e, rows_e, &tm);
        b_coo.push_block(rows_e, &cols_p, &gm);
        b_coo.push_block(&cols_p, rows_e, &gt);
        b_coo.push_block(&cols_p, &cols_p, &pp);
        if cell_cpl || cell_mu_cpl {
            has_coupling = true;
            if cell_mu_cpl {
                a1_coo.push_block(rows_e, rows_e, &a1_ee);
            }
            let neg: Vec<C> = mg.iter().map(|v| -*v).collect();
            a1_coo.push_block(rows_e, &cols_p, &neg);
            a1_coo.push_block(&cols_p, rows_e, &mgt);
        }
    }

    add_ibc_2d_waveguide(mesh, order, freq, &problem.pec, &problem.surface_impedance, &mut a_coo, &mut b_coo)?;
    let sheet_cpl = add_sheet2d_2d_waveguide(
        mesh,
        order,
        freq,
        &problem.pec,
        &problem.sheet_2d,
        us.slice_potential.as_deref(),
        &mut a_coo,
        &mut b_coo,
        &mut a1_coo,
    )?;
    has_coupling |= sheet_cpl;

    let a_full = a_coo.to_csr();
    let mut b_full = b_coo.to_csr();
    let mut non_pec_labels: Vec<i32> = problem.pmc.clone();
    for sb in &problem.sheet_2d {
        non_pec_labels.extend(sb.labels.iter().copied());
    }
    let (be, bn) = if problem.pec.is_empty()
        && problem.surface_impedance.is_empty()
        && non_pec_labels.is_empty()
        && problem.symmetry == Symmetry::None
    {
        let exterior = mesh.exterior_facets();
        let mut be: Vec<usize> = (0..mesh.n_facets())
            .filter(|&f| exterior[f])
            .map(|f| mesh.edge_index(mesh.facet(f)[0], mesh.facet(f)[1]))
            .collect();
        be.sort_unstable();
        be.dedup();
        let mut bn: Vec<usize> = (0..mesh.n_facets())
            .filter(|&f| exterior[f])
            .flat_map(|f| mesh.facet(f).iter().map(|&v| v as usize).collect::<Vec<_>>())
            .collect();
        bn.sort_unstable();
        bn.dedup();
        (be, bn)
    } else if problem.symmetry == Symmetry::None {
        pec_edges_and_nodes_2d(mesh, &problem.pec, &problem.surface_impedance, &non_pec_labels)
    } else {
        pec_edges_and_nodes_2d_with_symmetry(
            mesh,
            &problem.pec,
            &problem.surface_impedance,
            &non_pec_labels,
            problem.symmetry,
        )?
    };

    let mut fixed: Vec<usize> = be.clone();
    if order == 2 {
        fixed.extend(be.iter().map(|&e| ne + e));
    }
    fixed.extend(bn.iter().map(|&v| net + v));
    if order == 2 {
        fixed.extend(be.iter().map(|&e| net + nv + e));
    }
    let free = free_dofs(n_tot, &fixed);
    let scale: Vec<f64> = (0..n_tot).map(|i| if i < net { 1.0 } else { 1.0 / k0 }).collect();
    for r in 0..b_full.nrows() {
        for p in b_full.indptr()[r]..b_full.indptr()[r + 1] {
            let s = scale[r] * scale[b_full.indices()[p] as usize];
            b_full.values_mut()[p] *= s;
        }
    }
    let af = a_full.submatrix(&free, &free);
    let bf = b_full.submatrix(&free, &free);
    let a1 = if has_coupling {
        let mut a1 = a1_coo.to_csr();
        for r in 0..a1.nrows() {
            for p in a1.indptr()[r]..a1.indptr()[r + 1] {
                let s = scale[r] * scale[a1.indices()[p] as usize];
                a1.values_mut()[p] *= s;
            }
        }
        Some(a1.submatrix(&free, &free))
    } else {
        None
    };

    Ok(ModeAssemblyC { a: af, b: bf, a1, free, n_tot, net, k0, is_non_symmetric })
}

/// Filters PML modes, selects passive forward/backward propagation roots (`Im(beta) <= 0`),
/// applies biorthogonal / cross-power normalization (`|int (e_tilde_m^* x h_m) . y_hat dS| = 1`),
/// Gram–Schmidt orthogonalizes any degenerate clusters, and fixes deterministic spatial phase.
#[allow(clippy::too_many_arguments)]
fn finalize_slice_modes(
    us: &EmeUniqueSlice,
    problem: &Problem,
    freq: f64,
    k0: f64,
    has_pml: bool,
    pml_cells: &[bool],
    eps_base_w: &[C3],
    mui_phys_wafer: &[C3],
    mut mc: ModesC,
    mc_bwd: Option<ModesC>,
) -> Result<SolvedSliceModes, String> {
    let mesh = &us.mesh_2d;
    let order = problem.order;
    let slice_ibc_lossless = problem.surface_impedance.iter().all(|b| b.r_s.abs() <= 1e-15 && b.material.is_empty());
    let slice_is_hermitian = !has_pml
        && us.bend_radius_m.is_none()
        && slice_ibc_lossless
        && is_slice_hermitian(eps_base_w, mui_phys_wafer, &problem.sheet_2d, freq);
    let passive_lossless = slice_is_hermitian || (!has_pml && !problem.has_loss() && us.bend_radius_m.is_none());

    // Filter PML-dominated spurious modes when PML is active
    if has_pml && mc.neff.len() > 1 {
        let lf = compute_pml_layer_fractions(mesh, order, pml_cells, eps_base_w, &mc, k0);
        let keep: Vec<usize> = (0..mc.neff.len()).filter(|&i| lf[i] < 0.35).take(us.max_modes).collect();
        let sel = if keep.is_empty() { (0..mc.neff.len().min(us.max_modes)).collect::<Vec<_>>() } else { keep };
        mc.neff = sel.iter().map(|&i| mc.neff[i]).collect();
        mc.et = sel.iter().map(|&i| mc.et[i].clone()).collect();
        mc.phi = sel.iter().map(|&i| mc.phi[i].clone()).collect();
        mc.et_left = sel.iter().map(|&i| mc.et_left[i].clone()).collect();
        mc.phi_left = sel.iter().map(|&i| mc.phi_left[i].clone()).collect();
        mc.residuals = sel.iter().map(|&i| mc.residuals[i]).collect();
    } else if mc.neff.len() > us.max_modes {
        mc.neff.truncate(us.max_modes);
        mc.et.truncate(us.max_modes);
        mc.phi.truncate(us.max_modes);
        mc.et_left.truncate(us.max_modes);
        mc.phi_left.truncate(us.max_modes);
        mc.residuals.truncate(us.max_modes);
    }

    // Sort modes deterministically by descending Re(neff) for propagating modes, then ascending |Im(neff)| for evanescent modes
    let m_count = mc.neff.len();
    let mut order_idx: Vec<usize> = (0..m_count).collect();
    order_idx.sort_by(|&a, &b| {
        let (na, nb) = (mc.neff[a], mc.neff[b]);
        let prop_a = na.re.abs() >= na.im.abs() && na.re.abs() > 1e-4;
        let prop_b = nb.re.abs() >= nb.im.abs() && nb.re.abs() > 1e-4;
        prop_b.cmp(&prop_a).then_with(|| {
            if prop_a && prop_b {
                nb.re.abs().total_cmp(&na.re.abs())
            } else {
                na.im.abs().total_cmp(&nb.im.abs())
            }
        })
    });
    mc.neff = order_idx.iter().map(|&i| mc.neff[i]).collect();
    mc.et = order_idx.iter().map(|&i| mc.et[i].clone()).collect();
    mc.phi = order_idx.iter().map(|&i| mc.phi[i].clone()).collect();
    mc.et_left = order_idx.iter().map(|&i| mc.et_left[i].clone()).collect();
    mc.phi_left = order_idx.iter().map(|&i| mc.phi_left[i].clone()).collect();
    mc.residuals = order_idx.iter().map(|&i| mc.residuals[i]).collect();

    // Enforce forward (+y) root & passive evanescent decay Im(beta) <= 0
    let mut neff_fwd = Vec::with_capacity(m_count);
    let mut beta_fwd = Vec::with_capacity(m_count);
    for m in 0..m_count {
        let mut nf = mc.neff[m];
        let is_evan = nf.re.abs() < 1e-5 * nf.norm().max(1e-30) || nf.im.abs() > 20.0 * nf.re.abs();
        if is_evan {
            if nf.im > 0.0 {
                nf = -nf;
            }
            if passive_lossless {
                nf.re = 0.0;
            }
        } else {
            // Check real forward Poynting flux P_y; flip sign if P_y < 0 (e.g. negative root selected)
            let p_raw = slice_pair_integral(
                mesh,
                order,
                freq,
                pml_cells,
                mui_phys_wafer,
                &mc.et[m],
                &mc.et[m],
                &mc.phi[m],
                nf * k0,
            )
            .re;
            if (p_raw < -1e-14 && !mc.is_non_symmetric) || (nf.re < 0.0 && p_raw.abs() <= 1e-14) {
                nf = -nf;
            }
            let real_tol = if slice_is_hermitian { 1e-5 } else { 1e-7 };
            let im_cap = real_tol * nf.norm().max(1.0);
            if (passive_lossless && nf.im.abs() < im_cap) || (nf.im > 0.0 && nf.im < im_cap) {
                nf.im = 0.0;
            } else if nf.im > 0.0 {
                // Passive media cannot amplify along +y; flip to the decaying square root `-nf` (`(-nf)^2 == nf^2`)
                nf = -nf;
            }
        }
        neff_fwd.push(nf);
        beta_fwd.push(nf * k0);
    }

    // Backward effective indices (in passive decay convention Im(beta_bwd) <= 0 for exp(-j beta_bwd L_c))
    let mut neff_bwd = neff_fwd.clone();
    if let Some(bwd) = mc_bwd {
        let real_tol = if slice_is_hermitian { 1e-5 } else { 1e-7 };
        let mut bwd_pos: Vec<C> = bwd
            .neff
            .iter()
            .map(|&z| {
                let mut zp = -z;
                if passive_lossless && zp.im.abs() < real_tol * zp.norm().max(1.0) {
                    zp.im = 0.0;
                } else if zp.im > 0.0 {
                    zp.im = -zp.im.abs().min(0.0);
                }
                zp
            })
            .collect();
        bwd_pos.sort_by(|a, b| b.re.abs().total_cmp(&a.re.abs()));
        for m in 0..m_count.min(bwd_pos.len()) {
            let mut nb = bwd_pos[m];
            if nb.im > 0.0 {
                nb.im = 0.0;
            }
            neff_bwd[m] = nb;
        }
    }
    let beta_bwd: Vec<C> = neff_bwd.iter().map(|&nf| nf * k0).collect();

    // When lossless and free of complex-conjugate evanescent pairs (`Im(beta^2) == 0`), use Hermitian left mode = right mode (`et_left = et`);
    // when non-symmetric, PML/lossy, or complex-conjugate evanescent, keep true/complex-symmetric left eigenvectors.
    let has_complex_beta2 = beta_fwd.iter().any(|&b| (b * b).im.abs() > 1e-4 * b.norm_sqr().max(1e-30));
    let hermitian_media = slice_is_hermitian && !has_complex_beta2;
    if hermitian_media {
        mc.et_left = mc.et.clone();
        mc.phi_left = mc.phi.clone();
    }

    // Degenerate-cluster Gram–Schmidt orthogonalization under the cross-power pairing
    for m in 0..m_count {
        for k in 0..m {
            let rel_diff = (beta_fwd[m] - beta_fwd[k]).norm() / beta_fwd[k].norm().max(k0);
            if rel_diff < 1e-5 {
                let g_kk = slice_pair_integral(
                    mesh,
                    order,
                    freq,
                    pml_cells,
                    mui_phys_wafer,
                    &mc.et_left[k],
                    &mc.et[k],
                    &mc.phi[k],
                    beta_fwd[k],
                );
                if g_kk.norm() > 1e-24 {
                    let g_km = slice_pair_integral(
                        mesh,
                        order,
                        freq,
                        pml_cells,
                        mui_phys_wafer,
                        &mc.et_left[k],
                        &mc.et[m],
                        &mc.phi[m],
                        beta_fwd[m],
                    );
                    let alpha = g_km / g_kk;
                    for d in 0..mc.et[m].len() {
                        let ek = mc.et[k][d];
                        mc.et[m][d] -= alpha * ek;
                    }
                    for d in 0..mc.phi[m].len() {
                        let pk = mc.phi[k][d];
                        mc.phi[m][d] -= alpha * pk;
                    }
                    if hermitian_media {
                        mc.et_left[m] = mc.et[m].clone();
                        mc.phi_left[m] = mc.phi[m].clone();
                    } else if !mc.is_non_symmetric {
                        mc.et_left[m] = mc.et[m].iter().map(|z| z.conj()).collect();
                        mc.phi_left[m] = mc.phi[m].iter().map(|z| z.conj()).collect();
                    }
                }
            }
        }

        // Normalize mode m so |int (e_tilde_m^* x h_m) . y_hat dS| == 1 and fix deterministic spatial phase
        let phase_rot = canonical_mode_phase(mesh, order, pml_cells, &mc.et[m]);
        let g_mm = slice_pair_integral(
            mesh,
            order,
            freq,
            pml_cells,
            mui_phys_wafer,
            &mc.et_left[m],
            &mc.et[m],
            &mc.phi[m],
            beta_fwd[m],
        );
        let mag = g_mm.norm().max(1e-30).sqrt();
        let scale_r = phase_rot / mag;
        let scale_l = if !hermitian_media && !mc.is_non_symmetric { phase_rot.conj() / mag } else { phase_rot / mag };
        for v in &mut mc.et[m] {
            *v *= scale_r;
        }
        for v in &mut mc.phi[m] {
            *v *= scale_r;
        }
        for v in &mut mc.et_left[m] {
            *v *= scale_l;
        }
        for v in &mut mc.phi_left[m] {
            *v *= scale_l;
        }
    }

    let mut d_self = Vec::with_capacity(m_count);
    let mut d_self_tilde = Vec::with_capacity(m_count);
    for m in 0..m_count {
        let dm = slice_pair_integral(
            mesh,
            order,
            freq,
            pml_cells,
            mui_phys_wafer,
            &mc.et_left[m],
            &mc.et[m],
            &mc.phi[m],
            beta_fwd[m],
        );
        let dmt = slice_dual_pair_integral(
            mesh,
            order,
            freq,
            pml_cells,
            mui_phys_wafer,
            &mc.et[m],
            &mc.et_left[m],
            &mc.phi_left[m],
            beta_fwd[m],
        );
        d_self.push(dm);
        d_self_tilde.push(dmt);
    }

    Ok(SolvedSliceModes {
        mesh_2d: mesh.clone(),
        order,
        freq_hz: freq,
        bend_radius_m: us.bend_radius_m,
        pml_cells: pml_cells.to_vec(),
        mui_phys_wafer: mui_phys_wafer.to_vec(),
        neff_fwd,
        neff_bwd,
        beta_fwd,
        beta_bwd,
        et: mc.et,
        phi: mc.phi,
        et_left: mc.et_left,
        phi_left: mc.phi_left,
        d_self,
        d_self_tilde,
        residuals: mc.residuals,
        restarts: mc.restarts,
        converged: mc.converged,
        is_non_symmetric: mc.is_non_symmetric,
    })
}

fn is_slice_hermitian(eps_w: &[C3], mui_w: &[C3], sheets: &[Sheet2dBc], freq: f64) -> bool {
    let mat_herm = |arr: &[C3]| -> bool {
        arr.iter().all(|m| {
            let mut ok = true;
            for i in 0..3 {
                for j in 0..3 {
                    if (m[i][j] - m[j][i].conj()).norm() > 1e-13 {
                        ok = false;
                    }
                }
            }
            ok
        })
    };
    mat_herm(eps_w) && mat_herm(mui_w) && sheets.iter().all(|s| !s.has_loss(freq))
}

fn compute_pml_layer_fractions(
    mesh: &SimplexMesh,
    order: usize,
    pml_cells: &[bool],
    eps_w: &[C3],
    mc: &ModesC,
    k0: f64,
) -> Vec<f64> {
    let n_modes = mc.neff.len();
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let quad = quad_rule(2, 4);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, order)).collect();
    let nloc_e = if order == 1 { 3 } else { 8 };
    let nl = mesh.nloc(order);
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut u_phys = vec![0.0; n_modes];
    let mut u_pml = vec![0.0; n_modes];

    for ci in 0..mesh.n_cells() {
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        let er = eps_w[ci][0][0].re.abs().max(eps_w[ci][2][2].re.abs()).max(1.0);
        mesh.cell_dofs(ci, order, &mut pd);
        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
            let (nv_, _) = &tab[q];
            let wq = quad.w[q] * area * er;
            for j in 0..n_modes {
                let jb = c(0.0, 1.0) * (mc.neff[j] * k0);
                let (mut ex, mut ez) = (c(0.0, 0.0), c(0.0, 0.0));
                for s in 0..nloc_e {
                    let cf = mc.et[j][cd.gdofs[ci][s]];
                    ex += cf * w[s][0];
                    ez += cf * w[s][1];
                }
                let mut phi_q = c(0.0, 0.0);
                for a in 0..nl {
                    phi_q += mc.phi[j][pd[a]] * nv_[a];
                }
                let ey = jb * phi_q;
                let e2 = ex.norm_sqr() + ey.norm_sqr() + ez.norm_sqr();
                if pml_cells[ci] {
                    u_pml[j] += wq * e2;
                } else {
                    u_phys[j] += wq * e2;
                }
            }
        }
    }
    (0..n_modes)
        .map(|j| {
            let tot = u_phys[j] + u_pml[j];
            if tot > 0.0 {
                u_pml[j] / tot
            } else {
                0.0
            }
        })
        .collect()
}

/// Evaluates `int_{Omega_phys} (e_left^* x h_right) . y_hat dS` on `mesh`.
#[allow(clippy::too_many_arguments)]
fn slice_pair_integral(
    mesh: &SimplexMesh,
    order: usize,
    freq: f64,
    pml_cells: &[bool],
    mui_wafer: &[C3],
    et_left: &[C],
    et_right: &[C],
    phi_right: &[C],
    beta_right: C,
) -> C {
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let quad = quad_rule(2, 6);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, order)).collect();
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let omega = (2.0 * PI * freq).max(1e-300);
    let j_beta = c(0.0, 1.0) * beta_right;
    let h_pref = c(0.0, -1.0 / (omega * MU0));
    let ps = PROP_SIGN;

    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut gn = vec![[0.0; 2]; nl];
    let mut total = c(0.0, 0.0);

    for ci in 0..mesh.n_cells() {
        if pml_cells[ci] {
            continue;
        }
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        mesh.cell_dofs(ci, order, &mut pd);
        let mui_loc = wafer_to_local_mui(&mui_wafer[ci], ps);
        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
            let (_, dn) = &tab[q];
            for a in 0..nl {
                gn[a][0] = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
                gn[a][1] = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
            }
            let wq = quad.w[q] * area;
            let (mut el0, mut el1) = (c(0.0, 0.0), c(0.0, 0.0));
            let (mut er0, mut er1, mut curl_2d) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            for s in 0..nloc_e {
                let dof = cd.gdofs[ci][s];
                let cl = et_left[dof];
                let cr = et_right[dof];
                el0 += cl * w[s][0];
                el1 += cl * w[s][1];
                er0 += cr * w[s][0];
                er1 += cr * w[s][1];
                curl_2d += cr * curl[s];
            }
            let (mut dp0, mut dp1) = (c(0.0, 0.0), c(0.0, 0.0));
            for a in 0..nl {
                let pr = phi_right[pd[a]];
                dp0 += pr * gn[a][0];
                dp1 += pr * gn[a][1];
            }
            let curl_e = [j_beta * (dp1 + er1), -j_beta * (dp0 + er0), curl_2d];
            let mce = c3_mul_vec(&mui_loc, curl_e);
            let h0 = mce[0] * h_pref;
            let h1 = mce[1] * h_pref;
            // In wafer axes: (E_left^* x H_right) . y_hat = Ez_left^* Hx_right - Ex_left^* Hz_right = el1.conj() * h0 - el0.conj() * h1
            total += wq * (el1.conj() * h0 - el0.conj() * h1);
        }
    }
    total
}

/// Evaluates `int_{Omega_phys} (e_right x h_tilde_left^*) . y_hat dS` on `mesh`.
#[allow(clippy::too_many_arguments)]
fn slice_dual_pair_integral(
    mesh: &SimplexMesh,
    order: usize,
    freq: f64,
    pml_cells: &[bool],
    mui_wafer: &[C3],
    et_right: &[C],
    et_tilde: &[C],
    phi_tilde: &[C],
    beta_tilde_conj_source: C,
) -> C {
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let quad = quad_rule(2, 6);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, order)).collect();
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let omega = (2.0 * PI * freq).max(1e-300);
    // Left mode obeys adjoint equation with (eps^H, mu^H, beta^*); when Hermitian, h_tilde == h(beta)
    let j_beta = c(0.0, 1.0) * beta_tilde_conj_source.conj();
    let h_pref = c(0.0, -1.0 / (omega * MU0));
    let ps = PROP_SIGN;

    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut gn = vec![[0.0; 2]; nl];
    let mut total = c(0.0, 0.0);

    for ci in 0..mesh.n_cells() {
        if pml_cells[ci] {
            continue;
        }
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        mesh.cell_dofs(ci, order, &mut pd);
        let mui_adj_w = conj_transpose_c3(&mui_wafer[ci]);
        let mui_loc = wafer_to_local_mui(&mui_adj_w, ps);
        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
            let (_, dn) = &tab[q];
            for a in 0..nl {
                gn[a][0] = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
                gn[a][1] = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
            }
            let wq = quad.w[q] * area;
            let (mut er0, mut er1) = (c(0.0, 0.0), c(0.0, 0.0));
            let (mut et0, mut et1, mut curl_2d) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            for s in 0..nloc_e {
                let dof = cd.gdofs[ci][s];
                let cr = et_right[dof];
                let ct = et_tilde[dof];
                er0 += cr * w[s][0];
                er1 += cr * w[s][1];
                et0 += ct * w[s][0];
                et1 += ct * w[s][1];
                curl_2d += ct * curl[s];
            }
            let (mut dp0, mut dp1) = (c(0.0, 0.0), c(0.0, 0.0));
            for a in 0..nl {
                let pt = phi_tilde[pd[a]];
                dp0 += pt * gn[a][0];
                dp1 += pt * gn[a][1];
            }
            let curl_e = [j_beta * (dp1 + et1), -j_beta * (dp0 + et0), curl_2d];
            let mce = c3_mul_vec(&mui_loc, curl_e);
            let h0 = mce[0] * h_pref;
            let h1 = mce[1] * h_pref;
            total += wq * (er1 * h0.conj() - er0 * h1.conj());
        }
    }
    total
}

#[inline]
fn conj_transpose_c3(m: &C3) -> C3 {
    [
        [m[0][0].conj(), m[1][0].conj(), m[2][0].conj()],
        [m[0][1].conj(), m[1][1].conj(), m[2][1].conj()],
        [m[0][2].conj(), m[1][2].conj(), m[2][2].conj()],
    ]
}

#[inline]
fn wafer_to_local_mui(mui_w: &C3, ps: f64) -> C3 {
    let axes = [0usize, 2, 1];
    let mut loc = [[c(0.0, 0.0); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            loc[i][j] = mui_w[axes[i]][axes[j]];
        }
    }
    [
        [loc[0][0], loc[0][1], loc[0][2] * ps],
        [loc[1][0], loc[1][1], loc[1][2] * ps],
        [loc[2][0] * ps, loc[2][1] * ps, loc[2][2]],
    ]
}

/// Deterministic unit-modulus phase factor `e^{-j theta}` aligning the dominant spatial moment of `et` to positive real.
fn canonical_mode_phase(mesh: &SimplexMesh, order: usize, pml_cells: &[bool], et: &[C]) -> C {
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let quad = quad_rule(2, 4);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let (lo, hi) = mesh.bbox();
    let xc = 0.5 * (lo[0] + hi[0]);
    let zc = 0.5 * (lo[1] + hi[1]);
    let span = (hi[0] - lo[0]).hypot(hi[1] - lo[1]).max(1e-12);

    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut m0 = c(0.0, 0.0);
    let mut m1 = c(0.0, 0.0);
    let mut m2 = c(0.0, 0.0);
    let mut peak = c(0.0, 0.0);
    let mut area_tot = 0.0_f64;

    for ci in 0..mesh.n_cells() {
        if pml_cells[ci] {
            continue;
        }
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        let vs = mesh.cell(ci);
        let (p0, p1, p2) = (mesh.point(vs[0] as usize), mesh.point(vs[1] as usize), mesh.point(vs[2] as usize));
        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
            let wq = quad.w[q] * area;
            area_tot += wq;
            let (mut ex, mut ez) = (c(0.0, 0.0), c(0.0, 0.0));
            for s in 0..nloc_e {
                let cf = et[cd.gdofs[ci][s]];
                ex += cf * w[s][0];
                ez += cf * w[s][1];
            }
            let proj = ex + ez * std::f64::consts::FRAC_1_PI;
            let u = (l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0] - xc) / span;
            let v = (l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1] - zc) / span;
            m0 += proj * wq;
            m1 += proj * (wq * (u + 0.271_828_182_8 * v));
            m2 += proj * (wq * (u * u + 0.414_213_562_3 * u * v));
            if proj.norm_sqr() > peak.norm_sqr() {
                peak = proj;
            }
        }
    }
    let scale = area_tot.sqrt().max(1e-18);
    let chosen = if m0.norm() > 1e-4 * scale {
        m0
    } else if m1.norm() > 1e-5 * scale {
        m1
    } else if m2.norm() > 1e-6 * scale {
        m2
    } else {
        peak
    };
    if chosen.norm() > 1e-30 {
        chosen.conj() / chosen.norm()
    } else {
        c(1.0, 0.0)
    }
}

fn collect_tree_interface_pairs(node: &EmeCascadeNode, out: &mut Vec<(usize, usize)>) -> (usize, usize) {
    match node {
        EmeCascadeNode::Segment { unique_slice, .. } => (*unique_slice, *unique_slice),
        EmeCascadeNode::Composite { children } => {
            let (first, mut prev_last) = collect_tree_interface_pairs(&children[0], out);
            for ch in &children[1..] {
                let (ch_first, ch_last) = collect_tree_interface_pairs(ch, out);
                out.push((prev_last, ch_first));
                prev_last = ch_last;
            }
            (first, prev_last)
        }
        EmeCascadeNode::Periodic { unit_cell, num_reps, .. } => {
            let (first, last) = collect_tree_interface_pairs(unit_cell, out);
            // Always cache the periodic wrap interface `last -> first` so post-solve `eme_periodicity_sweep`
            // can evaluate arbitrary `N_rep >= 2` even when the nominal solve had `num_reps == 1`!
            let _ = num_reps;
            out.push((last, first));
            (first, last)
        }
    }
}

fn tree_end_slices(node: &EmeCascadeNode) -> (usize, usize) {
    match node {
        EmeCascadeNode::Segment { unique_slice, .. } => (*unique_slice, *unique_slice),
        EmeCascadeNode::Composite { children } => {
            let (first, _) = tree_end_slices(&children[0]);
            let (_, last) = tree_end_slices(&children[children.len() - 1]);
            (first, last)
        }
        EmeCascadeNode::Periodic { unit_cell, .. } => tree_end_slices(unit_cell),
    }
}

/// Computes the exact cross-power overlap matrix `O_{L,R}` and interface scattering matrix `S_if`
/// between slice `left_id` (`sm_l`) and slice `right_id` (`sm_r`).
fn compute_interface_overlap(
    left_id: usize,
    sm_l: &SolvedSliceModes,
    right_id: usize,
    sm_r: &SolvedSliceModes,
) -> Result<EmeInterfaceOverlap, String> {
    let m_l = sm_l.neff_fwd.len();
    let m_r = sm_r.neff_fwd.len();
    let d_left: Vec<[f64; 2]> = sm_l.d_self.iter().map(|z| [z.re, z.im]).collect();
    let d_right_tilde: Vec<[f64; 2]> = sm_r.d_self_tilde.iter().map(|z| [z.re, z.im]).collect();

    // Fast exact identity when both sides reference the same deduplicated slice!
    if left_id == right_id && m_l == m_r {
        let mut d_unit = vec![[1.0, 0.0]; m_l];
        for m in 0..m_l {
            let dm = sm_l.d_self[m];
            if (dm.re - 1.0).abs() >= 1e-3 || dm.im.abs() >= 1e-3 {
                d_unit[m] = [dm.re, dm.im];
            }
        }
        let mut v_lr = vec![vec![[0.0, 0.0]; m_r]; m_l];
        let mut w_lr = vec![vec![[0.0, 0.0]; m_r]; m_l];
        let mut o_cross = vec![vec![[0.0, 0.0]; m_r]; m_l];
        let mut s_if = vec![vec![c(0.0, 0.0); 2 * m_l]; 2 * m_l];
        for m in 0..m_l {
            v_lr[m][m] = d_unit[m];
            w_lr[m][m] = d_unit[m];
            o_cross[m][m] = [1.0, 0.0];
            s_if[m][m_l + m] = c(1.0, 0.0);
            s_if[m_l + m][m] = c(1.0, 0.0);
        }
        return Ok(EmeInterfaceOverlap {
            left_id,
            right_id,
            n_left: m_l,
            n_right: m_r,
            d_left: d_unit.clone(),
            d_right_tilde: d_unit,
            v_lr,
            w_lr,
            cross_power_overlap: o_cross,
            s_interface: c_mat_to_pairs(&s_if),
        });
    }

    let (v_mat, w_mat, o_mat) = if slices_geometrically_compatible(&sm_l.mesh_2d, &sm_r.mesh_2d) {
        overlap_shared_mesh(sm_l, sm_r)
    } else {
        overlap_nonmatching_meshes(sm_l, sm_r, false)?
    };

    let v_lr = c_mat_to_pairs(&v_mat);
    let w_lr = c_mat_to_pairs(&w_mat);
    let cross_power_overlap = c_mat_to_pairs(&o_mat);
    let s_if = build_interface_s_matrix(&d_left, &d_right_tilde, &v_lr, &w_lr, m_l, m_r)?;

    Ok(EmeInterfaceOverlap {
        left_id,
        right_id,
        n_left: m_l,
        n_right: m_r,
        d_left,
        d_right_tilde,
        v_lr,
        w_lr,
        cross_power_overlap,
        s_interface: c_mat_to_pairs(&s_if),
    })
}

/// Evaluates `(E, H, E_tilde, H_tilde)` at local barycentric coordinate `lam` of cell `ci` on a solved slice.
struct SliceCellEval<'a> {
    sm: &'a SolvedSliceModes,
    geo: femcore::fem::Geometry,
    cd: femcore::edge::CellData,
    omega: f64,
}

impl<'a> SliceCellEval<'a> {
    fn new(sm: &'a SolvedSliceModes) -> Self {
        let geo = geometry(&sm.mesh_2d);
        let cd = cell_data(&sm.mesh_2d, sm.order);
        let omega = (2.0 * PI * sm.freq_hz).max(1e-300);
        Self { sm, geo, cd, omega }
    }

    /// Evaluates `[Ex, Ez, Hx, Hz, Ex_tilde, Ez_tilde, Hx_tilde, Hz_tilde]` for every mode `0..M` at `(ci, lam)`.
    fn eval_all_modes(&self, ci: usize, lam: &[f64; 3], out: &mut Vec<[C; 8]>) {
        let m_count = self.sm.neff_fwd.len();
        out.clear();
        let order = self.sm.order;
        let nloc_e = if order == 1 { 3 } else { 8 };
        let nl = self.sm.mesh_2d.nloc(order);
        let g = tri_g(&self.geo, ci);
        let mut w = [[0.0; 2]; 8];
        let mut curl = [0.0; 8];
        edge_basis(lam, &g, &self.cd.sign[ci], order, &mut w[..nloc_e], &mut curl[..nloc_e]);
        let (_, dn) = shape(lam, order);
        let mut gn = [[0.0; 2]; 6];
        for a in 0..nl {
            gn[a][0] = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
            gn[a][1] = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
        }
        let mut pd = Vec::with_capacity(nl);
        self.sm.mesh_2d.cell_dofs(ci, order, &mut pd);

        let ps = PROP_SIGN;
        let h_pref = c(0.0, -1.0 / (self.omega * MU0));
        let mui_loc = wafer_to_local_mui(&self.sm.mui_phys_wafer[ci], ps);
        let mui_adj_loc = wafer_to_local_mui(&conj_transpose_c3(&self.sm.mui_phys_wafer[ci]), ps);

        for m in 0..m_count {
            let (mut ex, mut ez, mut curl_r) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            let (mut elx, mut elz, mut curl_l) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            for s in 0..nloc_e {
                let dof = self.cd.gdofs[ci][s];
                let cr = self.sm.et[m][dof];
                let cl = self.sm.et_left[m][dof];
                ex += cr * w[s][0];
                ez += cr * w[s][1];
                curl_r += cr * curl[s];
                elx += cl * w[s][0];
                elz += cl * w[s][1];
                curl_l += cl * curl[s];
            }
            let (mut dpr0, mut dpr1) = (c(0.0, 0.0), c(0.0, 0.0));
            let (mut dpl0, mut dpl1) = (c(0.0, 0.0), c(0.0, 0.0));
            for a in 0..nl {
                let pr = self.sm.phi[m][pd[a]];
                let pl = self.sm.phi_left[m][pd[a]];
                dpr0 += pr * gn[a][0];
                dpr1 += pr * gn[a][1];
                dpl0 += pl * gn[a][0];
                dpl1 += pl * gn[a][1];
            }
            let jb_r = c(0.0, 1.0) * self.sm.beta_fwd[m];
            let hr = c3_mul_vec(&mui_loc, [jb_r * (dpr1 + ez), -jb_r * (dpr0 + ex), curl_r]);
            let (hx, hz) = (hr[0] * h_pref, hr[1] * h_pref);

            let jb_l = c(0.0, 1.0) * self.sm.beta_fwd[m].conj();
            let hl = c3_mul_vec(&mui_adj_loc, [jb_l * (dpl1 + elz), -jb_l * (dpl0 + elx), curl_l]);
            let (hlx, hlz) = (hl[0] * h_pref, hl[1] * h_pref);

            out.push([ex, ez, hx, hz, elx, elz, hlx, hlz]);
        }
    }
}

/// Exact Gauss-Legendre interface overlap on two slices sharing an identical 2D triangulation.
fn overlap_shared_mesh(sm_l: &SolvedSliceModes, sm_r: &SolvedSliceModes) -> (Vec<Vec<C>>, Vec<Vec<C>>, Vec<Vec<C>>) {
    let m_l = sm_l.neff_fwd.len();
    let m_r = sm_r.neff_fwd.len();
    let eval_l = SliceCellEval::new(sm_l);
    let eval_r = SliceCellEval::new(sm_r);
    let quad = quad_rule(2, 6);
    let mut v_mat = vec![vec![c(0.0, 0.0); m_r]; m_l];
    let mut w_mat = vec![vec![c(0.0, 0.0); m_r]; m_l];
    let mut o_mat = vec![vec![c(0.0, 0.0); m_r]; m_l];
    let mut fl = Vec::with_capacity(m_l);
    let mut fr = Vec::with_capacity(m_r);

    for ci in 0..sm_l.mesh_2d.n_cells() {
        if sm_l.pml_cells[ci] || sm_r.pml_cells[ci] {
            continue;
        }
        let area = eval_l.geo.vol[ci];
        for (q, l) in quad.lam.iter().enumerate() {
            let lam = [l[0], l[1], l[2]];
            eval_l.eval_all_modes(ci, &lam, &mut fl);
            eval_r.eval_all_modes(ci, &lam, &mut fr);
            let wq = quad.w[q] * area;
            accumulate_overlap_point(&fl, &fr, wq, &mut v_mat, &mut w_mat, &mut o_mat);
        }
    }
    (v_mat, w_mat, o_mat)
}

/// Symmetric cross-mesh interface overlap between two arbitrary 2D slice triangulations via `femcore::transfer::TransferIndex`.
/// When `left_patch_only` is true (e.g. a sub-face `waveport` patch on `L`), integrates strictly over `sm_l` into `sm_r`.
fn overlap_nonmatching_meshes(
    sm_l: &SolvedSliceModes,
    sm_r: &SolvedSliceModes,
    left_patch_only: bool,
) -> Result<(Vec<Vec<C>>, Vec<Vec<C>>, Vec<Vec<C>>), String> {
    let m_l = sm_l.neff_fwd.len();
    let m_r = sm_r.neff_fwd.len();
    let eval_l = SliceCellEval::new(sm_l);
    let eval_r = SliceCellEval::new(sm_r);
    let empty_names = BTreeMap::new();
    let opts = TransferOptions { snap_rel: 0.0, parallel: false, ..TransferOptions::default() };
    let idx_r =
        TransferIndex::new(SourceMesh::Simplex(&sm_r.mesh_2d), &empty_names, opts).map_err(|e| e.to_string())?;
    let (fine_lam, fine_w) = fine_rule();

    let mut v_mat = vec![vec![c(0.0, 0.0); m_r]; m_l];
    let mut w_mat = vec![vec![c(0.0, 0.0); m_r]; m_l];
    let mut o_mat = vec![vec![c(0.0, 0.0); m_r]; m_l];
    let mut fl = Vec::with_capacity(m_l);
    let mut fr = Vec::with_capacity(m_r);

    let weight_l = if left_patch_only { 1.0 } else { 0.5 };
    for ci_l in 0..sm_l.mesh_2d.n_cells() {
        if sm_l.pml_cells[ci_l] {
            continue;
        }
        let area = eval_l.geo.vol[ci_l];
        let vs = sm_l.mesh_2d.cell(ci_l);
        let (p0, p1, p2) = (
            sm_l.mesh_2d.point(vs[0] as usize),
            sm_l.mesh_2d.point(vs[1] as usize),
            sm_l.mesh_2d.point(vs[2] as usize),
        );
        for (q, l) in fine_lam.iter().enumerate() {
            let xq = [l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0], l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1]];
            let Some(loc_r) = idx_r.locate(&xq, None) else { continue };
            if sm_r.pml_cells[loc_r.cell] {
                continue;
            }
            let lam_r = [loc_r.coords[0], loc_r.coords[1], loc_r.coords[2]];
            eval_l.eval_all_modes(ci_l, l, &mut fl);
            eval_r.eval_all_modes(loc_r.cell, &lam_r, &mut fr);
            accumulate_overlap_point(&fl, &fr, weight_l * fine_w[q] * area, &mut v_mat, &mut w_mat, &mut o_mat);
        }
    }

    if !left_patch_only {
        let idx_l =
            TransferIndex::new(SourceMesh::Simplex(&sm_l.mesh_2d), &empty_names, opts).map_err(|e| e.to_string())?;
        for ci_r in 0..sm_r.mesh_2d.n_cells() {
            if sm_r.pml_cells[ci_r] {
                continue;
            }
            let area = eval_r.geo.vol[ci_r];
            let vs = sm_r.mesh_2d.cell(ci_r);
            let (p0, p1, p2) = (
                sm_r.mesh_2d.point(vs[0] as usize),
                sm_r.mesh_2d.point(vs[1] as usize),
                sm_r.mesh_2d.point(vs[2] as usize),
            );
            for (q, l) in fine_lam.iter().enumerate() {
                let xq = [l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0], l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1]];
                let Some(loc_l) = idx_l.locate(&xq, None) else { continue };
                if sm_l.pml_cells[loc_l.cell] {
                    continue;
                }
                let lam_l = [loc_l.coords[0], loc_l.coords[1], loc_l.coords[2]];
                eval_l.eval_all_modes(loc_l.cell, &lam_l, &mut fl);
                eval_r.eval_all_modes(ci_r, l, &mut fr);
                accumulate_overlap_point(&fl, &fr, 0.5 * fine_w[q] * area, &mut v_mat, &mut w_mat, &mut o_mat);
            }
        }
    }

    Ok((v_mat, w_mat, o_mat))
}

#[inline]
fn accumulate_overlap_point(
    fl: &[[C; 8]],
    fr: &[[C; 8]],
    wq: f64,
    v_mat: &mut [Vec<C>],
    w_mat: &mut [Vec<C>],
    o_mat: &mut [Vec<C>],
) {
    let w_half = 0.5 * wq;
    for (m, lm) in fl.iter().enumerate() {
        let (ex_l, ez_l, hx_l, hz_l, elx_l, elz_l) = (lm[0], lm[1], lm[2], lm[3], lm[4], lm[5]);
        for (n, rn) in fr.iter().enumerate() {
            let (ex_r, ez_r, hx_r, hz_r, hlx_r, hlz_r) = (rn[0], rn[1], rn[2], rn[3], rn[6], rn[7]);
            // V_{LR}(m, n) = int (e_tilde_m^(L)* x h_n^(R)) . y_hat dS
            v_mat[m][n] += wq * (elz_l.conj() * hx_r - elx_l.conj() * hz_r);
            // W_{LR}(m, n) = int (e_m^(L) x h_tilde_n^(R)*) . y_hat dS
            w_mat[m][n] += wq * (ez_l * hlx_r.conj() - ex_l * hlz_r.conj());
            // O_{c,c+1}(m, n) = 0.5 int (E_m^(c) x H_n^(c+1)* + E_n^(c+1)* x H_m^(c)) . y_hat dS
            let e_l_x_h_r_conj = ez_l * hx_r.conj() - ex_l * hz_r.conj();
            let e_r_conj_x_h_l = ez_r.conj() * hx_l - ex_r.conj() * hz_l;
            o_mat[m][n] += w_half * (e_l_x_h_r_conj + e_r_conj_x_h_l);
        }
    }
}

/// Builds the `(k_l + k_r) x (k_l + k_r)` interface scattering matrix `S_if = M_-^{-1} M_+`
/// from truncated `(k_l, k_r)` sub-blocks of `(d_left, d_right_tilde, v_lr, w_lr)`:
/// `M_- = [D_L, V_LR; -W_LR^T, D_R_tilde]`, `M_+ = [D_L, V_LR; W_LR^T, -D_R_tilde]`.
fn build_interface_s_matrix(
    d_left: &[[f64; 2]],
    d_right_tilde: &[[f64; 2]],
    v_lr: &[Vec<[f64; 2]>],
    w_lr: &[Vec<[f64; 2]>],
    k_l: usize,
    k_r: usize,
) -> Result<Vec<Vec<C>>, String> {
    let n = k_l + k_r;
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut m_minus = vec![c(0.0, 0.0); n * n];
    let mut rhs_cols = vec![vec![c(0.0, 0.0); n]; n];

    for i in 0..k_l {
        let dl = c(d_left[i][0], d_left[i][1]);
        let dl_unit_pos = (dl.re - 1.0).abs() < 1e-4 && dl.im.abs() < 1e-4;
        m_minus[i * n + i] = if dl_unit_pos { c(1.0, 0.0) } else { dl };
        rhs_cols[i][i] = m_minus[i * n + i];
        for j in 0..k_r {
            let dr = c(d_right_tilde[j][0], d_right_tilde[j][1]);
            let dr_unit_pos = (dr.re - 1.0).abs() < 1e-4 && dr.im.abs() < 1e-4;
            let raw_v = c(v_lr[i][j][0], v_lr[i][j][1]);
            let raw_w = c(w_lr[i][j][0], w_lr[i][j][1]);
            // On lossless propagating pairs (D_ii = D_tilde_jj = +1) where |V_ij - W_ij^*| is pure quadrature noise,
            // enforce exact Hermitian adjointness W_ij = V_ij^* so S_if^H S_if == I to machine precision!
            let (vij, wij) =
                if dl_unit_pos && dr_unit_pos && (raw_v - raw_w.conj()).norm() < 1e-3 * (raw_v.norm() + 1.0) {
                    let u = c(0.5, 0.0) * (raw_v + raw_w.conj());
                    (u, u.conj())
                } else {
                    (raw_v, raw_w)
                };
            m_minus[i * n + (k_l + j)] = vij;
            rhs_cols[k_l + j][i] = vij;
            let row = k_l + j;
            m_minus[row * n + i] = -wij;
            rhs_cols[i][row] = wij;
        }
    }
    for j in 0..k_r {
        let row = k_l + j;
        let dr = c(d_right_tilde[j][0], d_right_tilde[j][1]);
        let dr_unit_pos = (dr.re - 1.0).abs() < 1e-4 && dr.im.abs() < 1e-4;
        let dr_eff = if dr_unit_pos { c(1.0, 0.0) } else { dr };
        m_minus[row * n + row] = dr_eff;
        rhs_cols[row][row] = -dr_eff;
    }

    let cols = solve_dense_complex_multi(&m_minus, &rhs_cols, n)
        .or_else(|| {
            let mut reg = m_minus.clone();
            for d in 0..n {
                reg[d * n + d] += c(1e-13, 0.0);
            }
            solve_dense_complex_multi(&reg, &rhs_cols, n)
        })
        .ok_or_else(|| "EME interface mode-matching matrix M_- is singular".to_string())?;

    let mut s_if = vec![vec![c(0.0, 0.0); n]; n];
    for col in 0..n {
        for row in 0..n {
            s_if[row][col] = cols[col][row];
        }
    }
    Ok(s_if)
}

/// Solves 2D modes on explicit `waveport` patches and builds port-to-slice interface overlaps.
fn build_explicit_waveport_couplings(
    asm: &EmeAssembly,
    first_sm: &SolvedSliceModes,
    first_u: usize,
    last_sm: &SolvedSliceModes,
    last_u: usize,
) -> Result<
    (
        Option<SolvedSliceModes>,
        Option<SolvedSliceModes>,
        Option<EmeInterfaceOverlap>,
        Option<EmeInterfaceOverlap>,
        Vec<(bool, usize)>,
        Vec<String>,
        Vec<f64>,
    ),
    String,
> {
    let mut in_patches: Vec<(usize, &EmeWavePortPatch)> = Vec::new();
    let mut out_patches: Vec<(usize, &EmeWavePortPatch)> = Vec::new();
    for (p_idx, patch) in asm.waveport_patches.iter().enumerate() {
        if patch.is_output {
            out_patches.push((p_idx, patch));
        } else {
            in_patches.push((p_idx, patch));
        }
    }

    let p_tot = asm.waveport_patches.len();
    let mut port_channels = vec![(false, 0usize); p_tot];
    let ports: Vec<String> = asm.waveport_patches.iter().map(|p| p.bc.name.clone()).collect();
    let port_z0_ohm: Vec<f64> = asm.waveport_patches.iter().map(|p| p.bc.z0_ohm.unwrap_or(asm.port_z0_ohm)).collect();

    // Solve each waveport patch's requested single mode, then compute overlap against slice 0 / slice C-1
    let solve_patch_single_mode = |patch: &EmeWavePortPatch| -> Result<SolvedSliceModes, String> {
        let mut sub_prob = asm.problem.clone();
        // Waveport patch rim uses label 1 = PEC and label 2 = PMC (matching `ports.rs`)
        sub_prob.pec = vec![1];
        sub_prob.pmc = vec![2];
        let m_req = patch.bc.mode.max(1);
        let (n_guess, k_patch) = if let Some(ng) = patch.bc.n_eff_guess {
            (ng, (m_req + 5).max(6))
        } else {
            let eps_w = sub_prob.materials.cell_full_complex_permittivity(&patch.mesh_2d, Some(1.0))?;
            let mur = sub_prob.materials.cell_values(&patch.mesh_2d, "mu_r", Some(1.0))?;
            let n_max = eps_w
                .iter()
                .zip(&mur)
                .map(|(e, &m)| {
                    let e_max = e[0][0].re.max(e[1][1].re).max(e[2][2].re).max(1.0);
                    (e_max * m.max(1.0)).sqrt()
                })
                .fold(1.0_f64, f64::max);
            (1.02 * n_max, (m_req + 3).max(4))
        };
        let us = EmeUniqueSlice {
            unique_id: 0,
            blake3_hash: String::new(),
            sample_y_m: if patch.is_output { asm.y_range_m[1] } else { asm.y_range_m[0] },
            mesh_2d: patch.mesh_2d.clone(),
            max_modes: k_patch,
            bend_radius_m: None,
            cell_quad4_temps: Vec::new(),
            cell_quad6_temps: Vec::new(),
            cell_avg_temps: Vec::new(),
            has_subcell_temp_grad: false,
            cell_quad4_e_dc: Vec::new(),
            cell_quad6_e_dc: Vec::new(),
            cell_avg_e_dc: Vec::new(),
            cell_quad4_strain: Vec::new(),
            cell_quad6_strain: Vec::new(),
            cell_avg_strain: Vec::new(),
            has_subcell_eo_pe_grad: false,
            slice_potential: None,
        };
        let mut sm = solve_unique_slice_modes(&us, &sub_prob, asm.freq_hz, n_guess, false)?;
        if sm.neff_fwd.len() < m_req {
            return Err(format!(
                "EME waveport '{}' requested mode {}, but only {} mode(s) converged",
                patch.bc.name,
                m_req,
                sm.neff_fwd.len()
            ));
        }
        let sel = if let Some(ng) = patch.bc.n_eff_guess {
            let mut order_by_guess: Vec<usize> = (0..sm.neff_fwd.len()).collect();
            order_by_guess.sort_by(|&a, &b| {
                let da = (sm.neff_fwd[a].re - ng).abs();
                let db = (sm.neff_fwd[b].re - ng).abs();
                da.total_cmp(&db)
            });
            order_by_guess[m_req - 1]
        } else {
            m_req - 1
        };
        sm.neff_fwd = vec![sm.neff_fwd[sel]];
        sm.neff_bwd = vec![sm.neff_bwd[sel]];
        sm.beta_fwd = vec![sm.beta_fwd[sel]];
        sm.beta_bwd = vec![sm.beta_bwd[sel]];
        sm.et = vec![sm.et[sel].clone()];
        sm.phi = vec![sm.phi[sel].clone()];
        sm.et_left = vec![sm.et_left[sel].clone()];
        sm.phi_left = vec![sm.phi_left[sel].clone()];
        sm.d_self = vec![sm.d_self[sel]];
        sm.d_self_tilde = vec![sm.d_self_tilde[sel]];
        sm.residuals = vec![sm.residuals[sel]];
        Ok(sm)
    };

    let input_ov = if !in_patches.is_empty() {
        let m_s = first_sm.neff_fwd.len();
        let n_in = in_patches.len();
        let mut d_left = Vec::with_capacity(n_in);
        let d_right_tilde: Vec<[f64; 2]> = first_sm.d_self_tilde.iter().map(|z| [z.re, z.im]).collect();
        let mut v_lr = Vec::with_capacity(n_in);
        let mut w_lr = Vec::with_capacity(n_in);
        let mut o_cross = Vec::with_capacity(n_in);
        for (ch, &(p_idx, patch)) in in_patches.iter().enumerate() {
            port_channels[p_idx] = (false, ch);
            let pm = solve_patch_single_mode(patch)?;
            d_left.push([pm.d_self[0].re, pm.d_self[0].im]);
            let (vm, wm, om) = overlap_nonmatching_meshes(&pm, first_sm, true)?;
            v_lr.push(vm[0].iter().map(|z| [z.re, z.im]).collect());
            w_lr.push(wm[0].iter().map(|z| [z.re, z.im]).collect());
            o_cross.push(om[0].iter().map(|z| [z.re, z.im]).collect());
        }
        let s_if = build_port_adapter_s_matrix(&o_cross, n_in, m_s, true);
        Some(EmeInterfaceOverlap {
            left_id: usize::MAX - 1,
            right_id: first_u,
            n_left: n_in,
            n_right: m_s,
            d_left,
            d_right_tilde,
            v_lr,
            w_lr,
            cross_power_overlap: o_cross,
            s_interface: c_mat_to_pairs(&s_if),
        })
    } else {
        None
    };

    let output_ov = if !out_patches.is_empty() {
        let m_s = last_sm.neff_fwd.len();
        let n_out = out_patches.len();
        let d_left: Vec<[f64; 2]> = last_sm.d_self.iter().map(|z| [z.re, z.im]).collect();
        let mut d_right_tilde = Vec::with_capacity(n_out);
        let mut v_lr = vec![vec![[0.0, 0.0]; n_out]; m_s];
        let mut w_lr = vec![vec![[0.0, 0.0]; n_out]; m_s];
        let mut o_cross = vec![vec![[0.0, 0.0]; n_out]; m_s];
        for (ch, &(p_idx, patch)) in out_patches.iter().enumerate() {
            port_channels[p_idx] = (true, ch);
            let pm = solve_patch_single_mode(patch)?;
            d_right_tilde.push([pm.d_self_tilde[0].re, pm.d_self_tilde[0].im]);
            let (vm_col, wm_col, om_col) = overlap_slice_to_patch(last_sm, &pm)?;
            for m in 0..m_s {
                v_lr[m][ch] = [vm_col[m].re, vm_col[m].im];
                w_lr[m][ch] = [wm_col[m].re, wm_col[m].im];
                o_cross[m][ch] = [om_col[m].re, om_col[m].im];
            }
        }
        let s_if = build_port_adapter_s_matrix(&o_cross, m_s, n_out, false);
        Some(EmeInterfaceOverlap {
            left_id: last_u,
            right_id: usize::MAX,
            n_left: m_s,
            n_right: n_out,
            d_left,
            d_right_tilde,
            v_lr,
            w_lr,
            cross_power_overlap: o_cross,
            s_interface: c_mat_to_pairs(&s_if),
        })
    } else {
        None
    };

    Ok((None, None, input_ov, output_ov, port_channels, ports, port_z0_ohm))
}

/// Builds a transparent matched modal port adapter scattering matrix (`S_{11} = 0, S_{22} = 0`)
/// projecting between `p` boundary `waveport` channels and `k` semi-infinite waveguide slice modes
/// without artificially reflecting unmonitored higher-order slice modes back into the device.
fn build_port_adapter_s_matrix(
    o_cross: &[Vec<[f64; 2]>],
    n_left: usize,
    n_right: usize,
    is_input: bool,
) -> Vec<Vec<C>> {
    let n = n_left + n_right;
    let mut s = vec![vec![c(0.0, 0.0); n]; n];
    if n_left == 0 || n_right == 0 {
        return s;
    }
    if is_input {
        let (p, k) = (n_left, n_right);
        let mut b = vec![vec![c(0.0, 0.0); p]; k];
        for ch in 0..p.min(o_cross.len()) {
            for m in 0..k.min(o_cross[ch].len()) {
                b[m][ch] = c(o_cross[ch][m][0], o_cross[ch][m][1]);
            }
        }
        clamp_tall_contractive(&mut b, k, p);
        for m in 0..k {
            for ch in 0..p {
                s[p + m][ch] = b[m][ch];
                s[ch][p + m] = b[m][ch].conj();
            }
        }
    } else {
        let (k, p) = (n_left, n_right);
        let mut b = vec![vec![c(0.0, 0.0); p]; k];
        for m in 0..k.min(o_cross.len()) {
            for ch in 0..p.min(o_cross[m].len()) {
                b[m][ch] = c(o_cross[m][ch][0], -o_cross[m][ch][1]);
            }
        }
        clamp_tall_contractive(&mut b, k, p);
        for m in 0..k {
            for ch in 0..p {
                s[k + ch][m] = b[m][ch].conj();
                s[m][k + ch] = b[m][ch];
            }
        }
    }
    s
}

/// Ensures `B` (`k x p`) is strictly contractive (`B^H B <= I_p`) by clamping any eigenvalue `> 1` of `G = B^H B`
/// via deterministic Hermitian Jacobi diagonalization (`B <- B V diag(min(1, lam^{-1/2})) V^H`).
fn clamp_tall_contractive(b: &mut [Vec<C>], k: usize, p: usize) {
    if k == 0 || p == 0 {
        return;
    }
    if p == 1 {
        let norm2: f64 = (0..k).map(|m| b[m][0].norm_sqr()).sum();
        if norm2 > 1.0 {
            let inv = 1.0 / norm2.sqrt();
            for row in b.iter_mut().take(k) {
                row[0] *= inv;
            }
        }
        return;
    }
    let mut g = vec![vec![c(0.0, 0.0); p]; p];
    for a in 0..p {
        for d in a..p {
            let mut acc = c(0.0, 0.0);
            for row in b.iter().take(k) {
                acc += row[a].conj() * row[d];
            }
            if a == d {
                g[a][a] = c(acc.re.max(0.0), 0.0);
            } else {
                g[a][d] = acc;
                g[d][a] = acc.conj();
            }
        }
    }
    let mut v = vec![vec![c(0.0, 0.0); p]; p];
    for i in 0..p {
        v[i][i] = c(1.0, 0.0);
    }
    for _sweep in 0..32 {
        let mut max_off = 0.0_f64;
        for pi in 0..p {
            for qi in (pi + 1)..p {
                let g_pq = g[pi][qi];
                let r = g_pq.norm();
                max_off = max_off.max(r);
                let scale = g[pi][pi].re.abs() + g[qi][qi].re.abs();
                if r <= 1e-15 * scale.max(1e-30) {
                    continue;
                }
                let phase = g_pq / r; // |phase| = 1, g_pq * conj(phase) = r > 0
                let dp = g[pi][pi].re;
                let dq = g[qi][qi].re;
                let tau = (dp - dq) / (2.0 * r);
                let t = if tau >= 0.0 {
                    1.0 / (tau + (1.0 + tau * tau).sqrt())
                } else {
                    -1.0 / (-tau + (1.0 + tau * tau).sqrt())
                };
                let cs = 1.0 / (1.0 + t * t).sqrt();
                let sn = t * cs;
                // Right-multiply U_step = [cs, -sn * phase; sn * conj(phase), cs] so U_step^H [dp, g_pq; conj(g_pq), dq] U_step is diagonal
                let u00 = c(cs, 0.0);
                let u01 = -phase * sn;
                let u10 = phase.conj() * sn;
                let u11 = c(cs, 0.0);
                for j in 0..p {
                    let (gp, gq) = (g[pi][j], g[qi][j]);
                    g[pi][j] = u00.conj() * gp + u10.conj() * gq;
                    g[qi][j] = u01.conj() * gp + u11.conj() * gq;
                }
                for i in 0..p {
                    let (gp, gq) = (g[i][pi], g[i][qi]);
                    g[i][pi] = gp * u00 + gq * u10;
                    g[i][qi] = gp * u01 + gq * u11;
                    let (vp, vq) = (v[i][pi], v[i][qi]);
                    v[i][pi] = vp * u00 + vq * u10;
                    v[i][qi] = vp * u01 + vq * u11;
                }
                g[pi][pi] = c(g[pi][pi].re.max(0.0), 0.0);
                g[qi][qi] = c(g[qi][qi].re.max(0.0), 0.0);
                g[pi][qi] = c(0.0, 0.0);
                g[qi][pi] = c(0.0, 0.0);
            }
        }
        if max_off <= 1e-15 {
            break;
        }
    }
    let scales: Vec<f64> = (0..p)
        .map(|i| {
            let lam = g[i][i].re;
            if lam > 1.0 {
                1.0 / lam.sqrt()
            } else {
                1.0
            }
        })
        .collect();
    if scales.iter().all(|&sc| (sc - 1.0).abs() == 0.0) {
        return;
    }
    let mut d_mat = vec![vec![c(0.0, 0.0); p]; p];
    for r in 0..p {
        for col in 0..p {
            let mut sum = c(0.0, 0.0);
            for e in 0..p {
                sum += v[r][e] * scales[e] * v[col][e].conj();
            }
            d_mat[r][col] = sum;
        }
    }
    for row in b.iter_mut().take(k) {
        let old = row.clone();
        for col in 0..p {
            let mut sum = c(0.0, 0.0);
            for q in 0..p {
                sum += old[q] * d_mat[q][col];
            }
            row[col] = sum;
        }
    }
}

/// Integrates slice `sm_slice` (`L`) against a single-mode waveport patch `sm_patch` (`R`) over `sm_patch`'s triangles.
fn overlap_slice_to_patch(
    sm_slice: &SolvedSliceModes,
    sm_patch: &SolvedSliceModes,
) -> Result<(Vec<C>, Vec<C>, Vec<C>), String> {
    let m_s = sm_slice.neff_fwd.len();
    let eval_l = SliceCellEval::new(sm_slice);
    let eval_r = SliceCellEval::new(sm_patch);
    let empty_names = BTreeMap::new();
    let opts = TransferOptions { snap_rel: 1e-4, parallel: false, ..TransferOptions::default() };
    let idx_l =
        TransferIndex::new(SourceMesh::Simplex(&sm_slice.mesh_2d), &empty_names, opts).map_err(|e| e.to_string())?;
    let (fine_lam, fine_w) = fine_rule();

    let mut v_mat = vec![vec![c(0.0, 0.0); 1]; m_s];
    let mut w_mat = vec![vec![c(0.0, 0.0); 1]; m_s];
    let mut o_mat = vec![vec![c(0.0, 0.0); 1]; m_s];
    let mut fl = Vec::with_capacity(m_s);
    let mut fr = Vec::with_capacity(1);

    for ci_r in 0..sm_patch.mesh_2d.n_cells() {
        let area = eval_r.geo.vol[ci_r];
        let vs = sm_patch.mesh_2d.cell(ci_r);
        let (p0, p1, p2) = (
            sm_patch.mesh_2d.point(vs[0] as usize),
            sm_patch.mesh_2d.point(vs[1] as usize),
            sm_patch.mesh_2d.point(vs[2] as usize),
        );
        for (q, l) in fine_lam.iter().enumerate() {
            let xq = [l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0], l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1]];
            let Some(loc_l) = idx_l.locate(&xq, None) else { continue };
            if sm_slice.pml_cells[loc_l.cell] {
                continue;
            }
            let lam_l = [loc_l.coords[0], loc_l.coords[1], loc_l.coords[2]];
            eval_l.eval_all_modes(loc_l.cell, &lam_l, &mut fl);
            eval_r.eval_all_modes(ci_r, l, &mut fr);
            accumulate_overlap_point(&fl, &fr, fine_w[q] * area, &mut v_mat, &mut w_mat, &mut o_mat);
        }
    }
    Ok((
        v_mat.into_iter().map(|r| r[0]).collect(),
        w_mat.into_iter().map(|r| r[0]).collect(),
        o_mat.into_iter().map(|r| r[0]).collect(),
    ))
}

// ------------------------------------------------------------------------------------------------------------- Redheffer Star Product & Sweeps

/// Scattering block `S` with `n_left` ports on the `-y` face and `n_right` ports on the `+y` face.
#[derive(Clone, Debug)]
struct ScatterBlock {
    n_left: usize,
    n_right: usize,
    /// `(n_left + n_right) x (n_left + n_right)` matrix partitioned as `[S11, S12; S21, S22]`.
    s: Vec<Vec<C>>,
}

/// Unconditionally stable propagation factor `exp(-j beta L)` with `Im(beta) <= 0` guaranteed (`|P| <= 1` for all `L >= 0`).
#[inline]
fn stable_prop_factor(beta: [f64; 2], length_m: f64) -> C {
    let l = length_m.max(0.0);
    // -j * (beta_re + j beta_im) * L = (beta_im * L) - j (beta_re * L); clamp beta_im <= 0 so real exponent <= 0!
    let decay = (beta[1].min(0.0) * l).max(-700.0).exp();
    let phase = -beta[0] * l;
    let (s, co) = phase.sin_cos();
    c(decay * co, decay * s)
}

fn propagation_block(beta_fwd: &[[f64; 2]], beta_bwd: &[[f64; 2]], k_modes: usize, length_m: f64) -> ScatterBlock {
    let m = k_modes.min(beta_fwd.len()).min(beta_bwd.len());
    let n = 2 * m;
    let mut s = vec![vec![c(0.0, 0.0); n]; n];
    for i in 0..m {
        // S21 = P_fwd (+y), S12 = P_bwd (-y)
        s[m + i][i] = stable_prop_factor(beta_fwd[i], length_m);
        s[i][m + i] = stable_prop_factor(beta_bwd[i], length_m);
    }
    ScatterBlock { n_left: m, n_right: m, s }
}

/// Redheffer star product `A * B` of `A` (`m1` left, `m2` right) and `B` (`m2` left, `m3` right) -> (`m1` left, `m3` right).
fn redheffer_star(a: &ScatterBlock, b: &ScatterBlock) -> Result<ScatterBlock, String> {
    if a.n_right != b.n_left {
        return Err(format!(
            "Redheffer star dimension mismatch: left block right-port count {} != right block left-port count {}",
            a.n_right, b.n_left
        ));
    }
    let (m1, m2, m3) = (a.n_left, a.n_right, b.n_right);
    let n_out = m1 + m3;
    if m2 == 0 {
        return Ok(ScatterBlock { n_left: m1, n_right: m3, s: vec![vec![c(0.0, 0.0); n_out]; n_out] });
    }

    // Extract S22^A (m2 x m2) and S11^B (m2 x m2), then form M1 = I - S22^A * S11^B and M2 = I - S11^B * S22^A
    let mut m1_mat = vec![c(0.0, 0.0); m2 * m2];
    let mut m2_mat = vec![c(0.0, 0.0); m2 * m2];
    for i in 0..m2 {
        m1_mat[i * m2 + i] = c(1.0, 0.0);
        m2_mat[i * m2 + i] = c(1.0, 0.0);
        for j in 0..m2 {
            let mut sum_ab = c(0.0, 0.0);
            let mut sum_ba = c(0.0, 0.0);
            for k in 0..m2 {
                sum_ab += a.s[m1 + i][m1 + k] * b.s[k][j];
                sum_ba += b.s[i][k] * a.s[m1 + k][m1 + j];
            }
            m1_mat[i * m2 + j] -= sum_ab;
            m2_mat[i * m2 + j] -= sum_ba;
        }
    }

    // Solve M2 * X1 = [S11^B * S21^A (m2 x m1) | S12^B (m2 x m3)] where M2 = I - S11^B * S22^A
    let mut rhs1 = vec![vec![c(0.0, 0.0); m2]; m1 + m3];
    for col in 0..m1 {
        for i in 0..m2 {
            let mut acc = c(0.0, 0.0);
            for k in 0..m2 {
                acc += b.s[i][k] * a.s[m1 + k][col];
            }
            rhs1[col][i] = acc;
        }
    }
    for col in 0..m3 {
        for i in 0..m2 {
            rhs1[m1 + col][i] = b.s[i][m2 + col];
        }
    }
    let sol1 = solve_with_reg(&m2_mat, &rhs1, m2)?;

    // Solve M1 * X2 = [S21^A (m2 x m1) | S22^A * S12^B (m2 x m3)] where M1 = I - S22^A * S11^B
    let mut rhs2 = vec![vec![c(0.0, 0.0); m2]; m1 + m3];
    for col in 0..m1 {
        for i in 0..m2 {
            rhs2[col][i] = a.s[m1 + i][col];
        }
    }
    for col in 0..m3 {
        for i in 0..m2 {
            let mut acc = c(0.0, 0.0);
            for k in 0..m2 {
                acc += a.s[m1 + i][m1 + k] * b.s[k][m2 + col];
            }
            rhs2[m1 + col][i] = acc;
        }
    }
    let sol2 = solve_with_reg(&m1_mat, &rhs2, m2)?;

    let mut s = vec![vec![c(0.0, 0.0); n_out]; n_out];
    // S11 = S11^A + S12^A * sol1[0..m1], S12 = S12^A * sol1[m1..m1+m3]
    for i in 0..m1 {
        for col in 0..m1 {
            let mut acc = a.s[i][col];
            for k in 0..m2 {
                acc += a.s[i][m1 + k] * sol1[col][k];
            }
            s[i][col] = acc;
        }
        for col in 0..m3 {
            let mut acc = c(0.0, 0.0);
            for k in 0..m2 {
                acc += a.s[i][m1 + k] * sol1[m1 + col][k];
            }
            s[i][m1 + col] = acc;
        }
    }
    // S21 = S21^B * sol2[0..m1], S22 = S22^B + S21^B * sol2[m1..m1+m3]
    for i in 0..m3 {
        for col in 0..m1 {
            let mut acc = c(0.0, 0.0);
            for k in 0..m2 {
                acc += b.s[m2 + i][k] * sol2[col][k];
            }
            s[m1 + i][col] = acc;
        }
        for col in 0..m3 {
            let mut acc = b.s[m2 + i][m2 + col];
            for k in 0..m2 {
                acc += b.s[m2 + i][k] * sol2[m1 + col][k];
            }
            s[m1 + i][m1 + col] = acc;
        }
    }

    Ok(ScatterBlock { n_left: m1, n_right: m3, s })
}

fn solve_with_reg(mat: &[C], rhs: &[Vec<C>], n: usize) -> Result<Vec<Vec<C>>, String> {
    if let Some(sol) = solve_dense_complex_multi(mat, rhs, n) {
        return Ok(sol);
    }
    let mut reg = mat.to_vec();
    for d in 0..n {
        reg[d * n + d] += c(1e-13, 0.0);
    }
    solve_dense_complex_multi(&reg, rhs, n)
        .ok_or_else(|| "Redheffer star internal denominator (I - S22 S11) is singular".into())
}

/// Binary exponentiation `block^{* exp}` (`exp >= 1`) in `O(log2 exp)` Redheffer star products.
fn redheffer_pow(block: &ScatterBlock, mut exp: usize) -> Result<ScatterBlock, String> {
    if exp == 0 || block.n_left != block.n_right {
        return Err("redheffer_pow requires exp >= 1 and a square left/right port count".into());
    }
    let mut base = block.clone();
    let mut acc: Option<ScatterBlock> = None;
    while exp > 0 {
        if (exp & 1) == 1 {
            acc = Some(match acc {
                None => base.clone(),
                Some(cur) => redheffer_star(&cur, &base)?,
            });
        }
        exp >>= 1;
        if exp > 0 {
            base = redheffer_star(&base, &base)?;
        }
    }
    acc.ok_or_else(|| "redheffer_pow internal state empty".into())
}

/// Evaluates the full `P x P` terminal scattering matrix from cached [`EmeData`] with optional overrides:
/// - `length_override`: `Some((segment_index, new_length_m))`
/// - `reps_override`: `Some((subgrid_index, new_num_reps))`
/// - `mode_cap`: `Some(max_modes_per_slice)`
fn evaluate_eme_terminal_s(
    data: &EmeData,
    length_override: Option<(usize, f64)>,
    reps_override: Option<(usize, usize)>,
    mode_cap: Option<usize>,
) -> Result<Vec<Vec<C>>, String> {
    let mut if_map: HashMap<(usize, usize, usize, usize), ScatterBlock> = HashMap::new();
    for ov in &data.interfaces {
        let k_l = mode_cap.map_or(ov.n_left, |mc| mc.min(ov.n_left)).max(1);
        let k_r = mode_cap.map_or(ov.n_right, |mc| mc.min(ov.n_right)).max(1);
        let s_mat = if k_l == ov.n_left && k_r == ov.n_right {
            pairs_to_c_mat(&ov.s_interface)
        } else {
            build_interface_s_matrix(&ov.d_left, &ov.d_right_tilde, &ov.v_lr, &ov.w_lr, k_l, k_r)?
        };
        if_map.insert((ov.left_id, ov.right_id, k_l, k_r), ScatterBlock { n_left: k_l, n_right: k_r, s: s_mat });
    }

    let (mut core_block, first_u, first_k, last_u, last_k) =
        eval_cascade_node(&data.tree, &data.unique_slices, &if_map, length_override, reps_override, mode_cap)?;

    if let Some(ref in_ov) = data.input_port_overlap {
        let s_in = if first_k == in_ov.n_right {
            pairs_to_c_mat(&in_ov.s_interface)
        } else {
            build_port_adapter_s_matrix(&in_ov.cross_power_overlap, in_ov.n_left, first_k, true)
        };
        let b_in = ScatterBlock { n_left: in_ov.n_left, n_right: first_k, s: s_in };
        core_block = redheffer_star(&b_in, &core_block)?;
    }
    let _ = (first_u, last_u);
    if let Some(ref out_ov) = data.output_port_overlap {
        let s_out = if last_k == out_ov.n_left {
            pairs_to_c_mat(&out_ov.s_interface)
        } else {
            build_port_adapter_s_matrix(&out_ov.cross_power_overlap, last_k, out_ov.n_right, false)
        };
        let b_out = ScatterBlock { n_left: last_k, n_right: out_ov.n_right, s: s_out };
        core_block = redheffer_star(&core_block, &b_out)?;
    }

    // Project onto `data.port_channels` (`P x P`)
    let p_count = data.port_channels.len();
    let mut s_term = vec![vec![c(0.0, 0.0); p_count]; p_count];
    for (r_p, &(r_out, r_ch)) in data.port_channels.iter().enumerate() {
        let r_idx = if r_out {
            if r_ch >= core_block.n_right {
                continue;
            }
            core_block.n_left + r_ch
        } else {
            if r_ch >= core_block.n_left {
                continue;
            }
            r_ch
        };
        for (c_p, &(c_out, c_ch)) in data.port_channels.iter().enumerate() {
            let c_idx = if c_out {
                if c_ch >= core_block.n_right {
                    continue;
                }
                core_block.n_left + c_ch
            } else {
                if c_ch >= core_block.n_left {
                    continue;
                }
                c_ch
            };
            s_term[r_p][c_p] = core_block.s[r_idx][c_idx];
        }
    }
    Ok(s_term)
}

#[allow(clippy::too_many_arguments)]
fn eval_cascade_node(
    node: &EmeCascadeNode,
    slices: &[EmeUniqueSliceSummary],
    if_map: &HashMap<(usize, usize, usize, usize), ScatterBlock>,
    length_override: Option<(usize, f64)>,
    reps_override: Option<(usize, usize)>,
    mode_cap: Option<usize>,
) -> Result<(ScatterBlock, usize, usize, usize, usize), String> {
    match node {
        EmeCascadeNode::Segment { segment_index, unique_slice, length_m, n_modes, .. } => {
            let us = &slices[*unique_slice];
            let k = mode_cap.map_or(*n_modes, |mc| mc.min(*n_modes)).min(us.n_modes).max(1);
            let len = match length_override {
                Some((idx, new_l)) if idx == *segment_index => new_l,
                _ => *length_m,
            };
            let blk = propagation_block(&us.beta_fwd, &us.beta_bwd, k, len);
            Ok((blk, *unique_slice, k, *unique_slice, k))
        }
        EmeCascadeNode::Composite { children } => {
            let (mut acc, first_u, first_k, mut prev_u, mut prev_k) =
                eval_cascade_node(&children[0], slices, if_map, length_override, reps_override, mode_cap)?;
            for ch in &children[1..] {
                let (ch_blk, ch_first_u, ch_first_k, ch_last_u, ch_last_k) =
                    eval_cascade_node(ch, slices, if_map, length_override, reps_override, mode_cap)?;
                let sif = if_map.get(&(prev_u, ch_first_u, prev_k, ch_first_k)).ok_or_else(|| {
                    format!("missing cached EME interface ({prev_u}, {ch_first_u}) for mode counts ({prev_k}, {ch_first_k})")
                })?;
                acc = redheffer_star(&acc, sif)?;
                acc = redheffer_star(&acc, &ch_blk)?;
                prev_u = ch_last_u;
                prev_k = ch_last_k;
            }
            Ok((acc, first_u, first_k, prev_u, prev_k))
        }
        EmeCascadeNode::Periodic { subgrid_index, unit_cell, num_reps } => {
            let reps = match reps_override {
                Some((idx, new_r)) if idx == *subgrid_index => new_r,
                _ => *num_reps,
            };
            if reps == 0 {
                return Err("EME PeriodicSubgrid repetition count must be >= 1".into());
            }
            let (cell_blk, first_u, first_k, last_u, last_k) =
                eval_cascade_node(unit_cell, slices, if_map, length_override, reps_override, mode_cap)?;
            if reps == 1 {
                return Ok((cell_blk, first_u, first_k, last_u, last_k));
            }
            let s_wrap = if_map
                .get(&(last_u, first_u, last_k, first_k))
                .ok_or_else(|| format!("missing cached EME periodic wrap interface ({last_u}, {first_u})"))?;
            let step = redheffer_star(&cell_blk, s_wrap)?;
            let pow_blk = redheffer_pow(&step, reps - 1)?;
            let full = redheffer_star(&pow_blk, &cell_blk)?;
            Ok((full, first_u, first_k, last_u, last_k))
        }
    }
}

/// Zero-re-solve analytical length sweep over segment `segment_index` (`lengths` in metres, `>= 0`).
pub fn eme_length_sweep(data: &EmeData, segment_index: usize, lengths: &[f64]) -> Result<EmeSweepReport, String> {
    let t0 = Instant::now();
    if segment_index >= data.total_segments {
        return Err(format!(
            "eme_length_sweep: segment_index {segment_index} out of bounds (cascade has {} segment(s))",
            data.total_segments
        ));
    }
    if lengths.is_empty() {
        return Err("eme_length_sweep requires at least one length value".into());
    }
    let mut points = Vec::with_capacity(lengths.len());
    for (i, &len) in lengths.iter().enumerate() {
        if !(len >= 0.0 && len.is_finite()) {
            return Err(format!("eme_length_sweep: lengths[{i}] = {len} must be finite and >= 0"));
        }
        let s_mat = evaluate_eme_terminal_s(data, Some((segment_index, len)), None, None)?;
        points.push(build_sweep_point(len, &s_mat));
    }
    Ok(EmeSweepReport {
        mode: "length".into(),
        target_index: segment_index,
        freq_hz: data.freq_hz,
        ports: data.ports.clone(),
        port_z0_ohm: data.port_z0_ohm.clone(),
        points,
        elapsed_s: t0.elapsed().as_secs_f64(),
    })
}

/// Zero-re-solve `O(log N_rep)` binary-exponentiation periodicity sweep over `PeriodicSubgrid` `subgrid_index`.
pub fn eme_periodicity_sweep(data: &EmeData, subgrid_index: usize, reps: &[usize]) -> Result<EmeSweepReport, String> {
    let t0 = Instant::now();
    if subgrid_index >= data.total_periodic_subgrids {
        return Err(format!(
            "eme_periodicity_sweep: subgrid_index {subgrid_index} out of bounds (cascade has {} PeriodicSubgrid node(s))",
            data.total_periodic_subgrids
        ));
    }
    if reps.is_empty() {
        return Err("eme_periodicity_sweep requires at least one repetition count".into());
    }
    let mut points = Vec::with_capacity(reps.len());
    for (i, &r) in reps.iter().enumerate() {
        if r == 0 {
            return Err(format!("eme_periodicity_sweep: reps[{i}] must be >= 1"));
        }
        let s_mat = evaluate_eme_terminal_s(data, None, Some((subgrid_index, r)), None)?;
        points.push(build_sweep_point(r as f64, &s_mat));
    }
    Ok(EmeSweepReport {
        mode: "periodicity".into(),
        target_index: subgrid_index,
        freq_hz: data.freq_hz,
        ports: data.ports.clone(),
        port_z0_ohm: data.port_z0_ohm.clone(),
        points,
        elapsed_s: t0.elapsed().as_secs_f64(),
    })
}

/// Zero-re-solve modal truncation convergence sweep truncating every slice to `1 <= k <= M_max` modes.
pub fn eme_mode_sweep(data: &EmeData, mode_counts: &[usize]) -> Result<EmeSweepReport, String> {
    let t0 = Instant::now();
    if mode_counts.is_empty() {
        return Err("eme_mode_sweep requires at least one mode count".into());
    }
    let mut points = Vec::with_capacity(mode_counts.len());
    for (i, &k) in mode_counts.iter().enumerate() {
        if k == 0 {
            return Err(format!("eme_mode_sweep: mode_counts[{i}] must be >= 1"));
        }
        let s_mat = evaluate_eme_terminal_s(data, None, None, Some(k))?;
        points.push(build_sweep_point(k as f64, &s_mat));
    }
    Ok(EmeSweepReport {
        mode: "modes".into(),
        target_index: 0,
        freq_hz: data.freq_hz,
        ports: data.ports.clone(),
        port_z0_ohm: data.port_z0_ohm.clone(),
        points,
        elapsed_s: t0.elapsed().as_secs_f64(),
    })
}

fn build_sweep_point(param: f64, s_mat: &[Vec<C>]) -> EmeSweepPoint {
    let (sigma_max, recip_err) = s_matrix_passivity_and_reciprocity(s_mat);
    let p = s_mat.len();
    let power_from_port0 = if p > 0 { (0..p).map(|r| s_mat[r][0].norm_sqr()).collect() } else { Vec::new() };
    EmeSweepPoint {
        parameter: param,
        s_matrix: c_mat_to_pairs(s_mat),
        power_from_port0,
        passivity_max_sigma: sigma_max,
        reciprocity_rel_err: recip_err,
    }
}

fn c_mat_to_pairs(m: &[Vec<C>]) -> Vec<Vec<[f64; 2]>> {
    m.iter().map(|r| r.iter().map(|z| [z.re, z.im]).collect()).collect()
}

fn pairs_to_c_mat(m: &[Vec<[f64; 2]>]) -> Vec<Vec<C>> {
    m.iter().map(|r| r.iter().map(|&[re, im]| c(re, im)).collect()).collect()
}

// ------------------------------------------------------------------------------------------------------------- Analytical Redheffer Sensitivities (FR-37 + FR-40)

fn zero_scatter_block(n_left: usize, n_right: usize) -> ScatterBlock {
    let n = n_left + n_right;
    ScatterBlock { n_left, n_right, s: vec![vec![c(0.0, 0.0); n]; n] }
}

fn propagation_block_and_deriv(
    beta_fwd: &[[f64; 2]],
    beta_bwd: &[[f64; 2]],
    k_modes: usize,
    length_m: f64,
    d_length: f64,
) -> (ScatterBlock, ScatterBlock) {
    let m = k_modes.min(beta_fwd.len()).min(beta_bwd.len());
    let n = 2 * m;
    let mut s = vec![vec![c(0.0, 0.0); n]; n];
    let mut ds = vec![vec![c(0.0, 0.0); n]; n];
    for i in 0..m {
        let pf = stable_prop_factor(beta_fwd[i], length_m);
        let pb = stable_prop_factor(beta_bwd[i], length_m);
        s[m + i][i] = pf;
        s[i][m + i] = pb;
        if d_length != 0.0 {
            // P(L) = exp((min(beta_im, 0) - j beta_re) L) => dP/dL = -j beta_passive P(L)
            let gf = c(beta_fwd[i][1].min(0.0), -beta_fwd[i][0]) * d_length;
            let gb = c(beta_bwd[i][1].min(0.0), -beta_bwd[i][0]) * d_length;
            ds[m + i][i] = gf * pf;
            ds[i][m + i] = gb * pb;
        }
    }
    (ScatterBlock { n_left: m, n_right: m, s }, ScatterBlock { n_left: m, n_right: m, s: ds })
}

/// Exact analytical tangent of [`redheffer_star`]: returns `(A * B, d(A * B))`.
fn redheffer_star_tangent(
    a: &ScatterBlock,
    da: &ScatterBlock,
    b: &ScatterBlock,
    db: &ScatterBlock,
) -> Result<(ScatterBlock, ScatterBlock), String> {
    if a.n_right != b.n_left {
        return Err(format!(
            "Redheffer star dimension mismatch: left block right-port count {} != right block left-port count {}",
            a.n_right, b.n_left
        ));
    }
    let (m1, m2, m3) = (a.n_left, a.n_right, b.n_right);
    let n_out = m1 + m3;
    if m2 == 0 {
        return Ok((zero_scatter_block(m1, m3), zero_scatter_block(m1, m3)));
    }

    let mut m1_mat = vec![c(0.0, 0.0); m2 * m2];
    let mut m2_mat = vec![c(0.0, 0.0); m2 * m2];
    for i in 0..m2 {
        m1_mat[i * m2 + i] = c(1.0, 0.0);
        m2_mat[i * m2 + i] = c(1.0, 0.0);
        for j in 0..m2 {
            let mut sum_ab = c(0.0, 0.0);
            let mut sum_ba = c(0.0, 0.0);
            for k in 0..m2 {
                sum_ab += a.s[m1 + i][m1 + k] * b.s[k][j];
                sum_ba += b.s[i][k] * a.s[m1 + k][m1 + j];
            }
            m1_mat[i * m2 + j] -= sum_ab;
            m2_mat[i * m2 + j] -= sum_ba;
        }
    }

    let mut rhs1 = vec![vec![c(0.0, 0.0); m2]; n_out];
    let mut rhs2 = vec![vec![c(0.0, 0.0); m2]; n_out];
    for col in 0..m1 {
        for i in 0..m2 {
            let mut acc = c(0.0, 0.0);
            for k in 0..m2 {
                acc += b.s[i][k] * a.s[m1 + k][col];
            }
            rhs1[col][i] = acc;
            rhs2[col][i] = a.s[m1 + i][col];
        }
    }
    for col in 0..m3 {
        for i in 0..m2 {
            rhs1[m1 + col][i] = b.s[i][m2 + col];
            let mut acc = c(0.0, 0.0);
            for k in 0..m2 {
                acc += a.s[m1 + i][m1 + k] * b.s[k][m2 + col];
            }
            rhs2[m1 + col][i] = acc;
        }
    }
    let v = solve_with_reg(&m2_mat, &rhs1, m2)?;
    let u = solve_with_reg(&m1_mat, &rhs2, m2)?;

    let mut s = vec![vec![c(0.0, 0.0); n_out]; n_out];
    for i in 0..m1 {
        for col in 0..m1 {
            let mut acc = a.s[i][col];
            for k in 0..m2 {
                acc += a.s[i][m1 + k] * v[col][k];
            }
            s[i][col] = acc;
        }
        for col in 0..m3 {
            let mut acc = c(0.0, 0.0);
            for k in 0..m2 {
                acc += a.s[i][m1 + k] * v[m1 + col][k];
            }
            s[i][m1 + col] = acc;
        }
    }
    for i in 0..m3 {
        for col in 0..m1 {
            let mut acc = c(0.0, 0.0);
            for k in 0..m2 {
                acc += b.s[m2 + i][k] * u[col][k];
            }
            s[m1 + i][col] = acc;
        }
        for col in 0..m3 {
            let mut acc = b.s[m2 + i][m2 + col];
            for k in 0..m2 {
                acc += b.s[m2 + i][k] * u[m1 + col][k];
            }
            s[m1 + i][m1 + col] = acc;
        }
    }

    // Fast skip when both incoming block derivatives are identically zero
    let da_zero = da.s.iter().flatten().all(|z| z.re == 0.0 && z.im == 0.0);
    let db_zero = db.s.iter().flatten().all(|z| z.re == 0.0 && z.im == 0.0);
    if da_zero && db_zero {
        return Ok((ScatterBlock { n_left: m1, n_right: m3, s }, zero_scatter_block(m1, m3)));
    }

    let mut delta_u = vec![vec![c(0.0, 0.0); m2]; n_out];
    let mut delta_v = vec![vec![c(0.0, 0.0); m2]; n_out];
    for col in 0..n_out {
        for i in 0..m2 {
            let mut du_acc = if col < m1 { da.s[m1 + i][col] } else { c(0.0, 0.0) };
            let mut dv_acc = if col >= m1 { db.s[i][m2 + (col - m1)] } else { c(0.0, 0.0) };
            for k in 0..m2 {
                du_acc += da.s[m1 + i][m1 + k] * v[col][k];
                dv_acc += db.s[i][k] * u[col][k];
            }
            delta_u[col][i] = du_acc;
            delta_v[col][i] = dv_acc;
        }
    }

    let mut drhs_v = vec![vec![c(0.0, 0.0); m2]; n_out];
    let mut drhs_u = vec![vec![c(0.0, 0.0); m2]; n_out];
    for col in 0..n_out {
        for i in 0..m2 {
            let mut bv = delta_v[col][i];
            let mut au = delta_u[col][i];
            for k in 0..m2 {
                bv += b.s[i][k] * delta_u[col][k];
                au += a.s[m1 + i][m1 + k] * delta_v[col][k];
            }
            drhs_v[col][i] = bv;
            drhs_u[col][i] = au;
        }
    }
    let dv = solve_with_reg(&m2_mat, &drhs_v, m2)?;
    let du = solve_with_reg(&m1_mat, &drhs_u, m2)?;

    let mut ds = vec![vec![c(0.0, 0.0); n_out]; n_out];
    for i in 0..m1 {
        for col in 0..m1 {
            let mut acc = da.s[i][col];
            for k in 0..m2 {
                acc += da.s[i][m1 + k] * v[col][k] + a.s[i][m1 + k] * dv[col][k];
            }
            ds[i][col] = acc;
        }
        for col in 0..m3 {
            let mut acc = c(0.0, 0.0);
            for k in 0..m2 {
                acc += da.s[i][m1 + k] * v[m1 + col][k] + a.s[i][m1 + k] * dv[m1 + col][k];
            }
            ds[i][m1 + col] = acc;
        }
    }
    for i in 0..m3 {
        for col in 0..m1 {
            let mut acc = c(0.0, 0.0);
            for k in 0..m2 {
                acc += db.s[m2 + i][k] * u[col][k] + b.s[m2 + i][k] * du[col][k];
            }
            ds[m1 + i][col] = acc;
        }
        for col in 0..m3 {
            let mut acc = db.s[m2 + i][m2 + col];
            for k in 0..m2 {
                acc += db.s[m2 + i][k] * u[m1 + col][k] + b.s[m2 + i][k] * du[m1 + col][k];
            }
            ds[m1 + i][m1 + col] = acc;
        }
    }

    Ok((ScatterBlock { n_left: m1, n_right: m3, s }, ScatterBlock { n_left: m1, n_right: m3, s: ds }))
}

fn redheffer_pow_tangent(
    block: &ScatterBlock,
    dblock: &ScatterBlock,
    mut exp: usize,
) -> Result<(ScatterBlock, ScatterBlock), String> {
    if exp == 0 || block.n_left != block.n_right {
        return Err("redheffer_pow_tangent requires exp >= 1 and a square left/right port count".into());
    }
    let mut base = (block.clone(), dblock.clone());
    let mut acc: Option<(ScatterBlock, ScatterBlock)> = None;
    while exp > 0 {
        if (exp & 1) == 1 {
            acc = Some(match acc {
                None => base.clone(),
                Some((cur_b, cur_db)) => redheffer_star_tangent(&cur_b, &cur_db, &base.0, &base.1)?,
            });
        }
        exp >>= 1;
        if exp > 0 {
            let (next_b, next_db) = redheffer_star_tangent(&base.0, &base.1, &base.0, &base.1)?;
            base = (next_b, next_db);
        }
    }
    acc.ok_or_else(|| "redheffer_pow_tangent internal state empty".into())
}

impl EmeData {
    /// Returns the nominal segment propagation lengths `L_c` (`0..self.total_segments`) stored in `self.tree`.
    pub fn segment_lengths(&self) -> Vec<f64> {
        let mut out = vec![0.0; self.total_segments];
        fn visit(node: &EmeCascadeNode, out: &mut [f64]) {
            match node {
                EmeCascadeNode::Segment { segment_index, length_m, .. } => {
                    if let Some(slot) = out.get_mut(*segment_index) {
                        *slot = *length_m;
                    }
                }
                EmeCascadeNode::Composite { children } => {
                    for ch in children {
                        visit(ch, out);
                    }
                }
                EmeCascadeNode::Periodic { unit_cell, .. } => visit(unit_cell, out),
            }
        }
        visit(&self.tree, &mut out);
        out
    }

    /// Evaluates both the `P x P` terminal scattering matrix `S` and its exact directional tangent
    /// `dS = sum_c (dS / dL_c) * d_lengths[c]` via `dP_c/dL_c = -j beta_c P_c` and `O(log N)` Redheffer tangents.
    pub(crate) fn evaluate_terminal_s_and_grad(
        &self,
        length_overrides: Option<&[f64]>,
        d_lengths: &[f64],
    ) -> Result<(Vec<Vec<C>>, Vec<Vec<C>>), String> {
        let mut if_map: HashMap<(usize, usize, usize, usize), ScatterBlock> = HashMap::new();
        for ov in &self.interfaces {
            let k_l = ov.n_left.max(1);
            let k_r = ov.n_right.max(1);
            if_map.insert(
                (ov.left_id, ov.right_id, k_l, k_r),
                ScatterBlock { n_left: k_l, n_right: k_r, s: pairs_to_c_mat(&ov.s_interface) },
            );
        }

        let (mut core_block, mut dcore_block, _first_u, first_k, _last_u, last_k) =
            eval_cascade_node_tangent(&self.tree, &self.unique_slices, &if_map, length_overrides, d_lengths)?;

        if let Some(ref in_ov) = self.input_port_overlap {
            let s_in = if first_k == in_ov.n_right {
                pairs_to_c_mat(&in_ov.s_interface)
            } else {
                build_port_adapter_s_matrix(&in_ov.cross_power_overlap, in_ov.n_left, first_k, true)
            };
            let b_in = ScatterBlock { n_left: in_ov.n_left, n_right: first_k, s: s_in };
            let db_in = zero_scatter_block(in_ov.n_left, first_k);
            let (nb, ndb) = redheffer_star_tangent(&b_in, &db_in, &core_block, &dcore_block)?;
            core_block = nb;
            dcore_block = ndb;
        }
        if let Some(ref out_ov) = self.output_port_overlap {
            let s_out = if last_k == out_ov.n_left {
                pairs_to_c_mat(&out_ov.s_interface)
            } else {
                build_port_adapter_s_matrix(&out_ov.cross_power_overlap, last_k, out_ov.n_right, false)
            };
            let b_out = ScatterBlock { n_left: last_k, n_right: out_ov.n_right, s: s_out };
            let db_out = zero_scatter_block(last_k, out_ov.n_right);
            let (nb, ndb) = redheffer_star_tangent(&core_block, &dcore_block, &b_out, &db_out)?;
            core_block = nb;
            dcore_block = ndb;
        }

        let p_count = self.port_channels.len();
        let mut s_term = vec![vec![c(0.0, 0.0); p_count]; p_count];
        let mut ds_term = vec![vec![c(0.0, 0.0); p_count]; p_count];
        for (r_p, &(r_out, r_ch)) in self.port_channels.iter().enumerate() {
            let r_idx = if r_out {
                if r_ch >= core_block.n_right {
                    continue;
                }
                core_block.n_left + r_ch
            } else {
                if r_ch >= core_block.n_left {
                    continue;
                }
                r_ch
            };
            for (c_p, &(c_out, c_ch)) in self.port_channels.iter().enumerate() {
                let c_idx = if c_out {
                    if c_ch >= core_block.n_right {
                        continue;
                    }
                    core_block.n_left + c_ch
                } else {
                    if c_ch >= core_block.n_left {
                        continue;
                    }
                    c_ch
                };
                s_term[r_p][c_p] = core_block.s[r_idx][c_idx];
                ds_term[r_p][c_p] = dcore_block.s[r_idx][c_idx];
            }
        }
        Ok((s_term, ds_term))
    }
}

fn eval_cascade_node_tangent(
    node: &EmeCascadeNode,
    slices: &[EmeUniqueSliceSummary],
    if_map: &HashMap<(usize, usize, usize, usize), ScatterBlock>,
    length_overrides: Option<&[f64]>,
    d_lengths: &[f64],
) -> Result<(ScatterBlock, ScatterBlock, usize, usize, usize, usize), String> {
    match node {
        EmeCascadeNode::Segment { segment_index, unique_slice, length_m, n_modes, .. } => {
            let us = &slices[*unique_slice];
            let k = (*n_modes).min(us.n_modes).max(1);
            let len = length_overrides.and_then(|ov| ov.get(*segment_index).copied()).unwrap_or(*length_m);
            let dlen = d_lengths.get(*segment_index).copied().unwrap_or(0.0);
            let (blk, dblk) = propagation_block_and_deriv(&us.beta_fwd, &us.beta_bwd, k, len, dlen);
            Ok((blk, dblk, *unique_slice, k, *unique_slice, k))
        }
        EmeCascadeNode::Composite { children } => {
            let (mut acc, mut dacc, first_u, first_k, mut prev_u, mut prev_k) =
                eval_cascade_node_tangent(&children[0], slices, if_map, length_overrides, d_lengths)?;
            for ch in &children[1..] {
                let (ch_blk, ch_dblk, ch_first_u, ch_first_k, ch_last_u, ch_last_k) =
                    eval_cascade_node_tangent(ch, slices, if_map, length_overrides, d_lengths)?;
                let sif = if_map.get(&(prev_u, ch_first_u, prev_k, ch_first_k)).ok_or_else(|| {
                    format!(
                        "missing cached EME interface ({prev_u}, {ch_first_u}) for mode counts ({prev_k}, {ch_first_k})"
                    )
                })?;
                let dsif = zero_scatter_block(sif.n_left, sif.n_right);
                let (a1, da1) = redheffer_star_tangent(&acc, &dacc, sif, &dsif)?;
                let (a2, da2) = redheffer_star_tangent(&a1, &da1, &ch_blk, &ch_dblk)?;
                acc = a2;
                dacc = da2;
                prev_u = ch_last_u;
                prev_k = ch_last_k;
            }
            Ok((acc, dacc, first_u, first_k, prev_u, prev_k))
        }
        EmeCascadeNode::Periodic { unit_cell, num_reps, .. } => {
            let reps = *num_reps;
            if reps == 0 {
                return Err("EME PeriodicSubgrid repetition count must be >= 1".into());
            }
            let (cell_blk, cell_dblk, first_u, first_k, last_u, last_k) =
                eval_cascade_node_tangent(unit_cell, slices, if_map, length_overrides, d_lengths)?;
            if reps == 1 {
                return Ok((cell_blk, cell_dblk, first_u, first_k, last_u, last_k));
            }
            let s_wrap = if_map
                .get(&(last_u, first_u, last_k, first_k))
                .ok_or_else(|| format!("missing cached EME periodic wrap interface ({last_u}, {first_u})"))?;
            let ds_wrap = zero_scatter_block(s_wrap.n_left, s_wrap.n_right);
            let (step, dstep) = redheffer_star_tangent(&cell_blk, &cell_dblk, s_wrap, &ds_wrap)?;
            let (pow_blk, dpow_blk) = redheffer_pow_tangent(&step, &dstep, reps - 1)?;
            let (full, dfull) = redheffer_star_tangent(&pow_blk, &dpow_blk, &cell_blk, &cell_dblk)?;
            Ok((full, dfull, first_u, first_k, last_u, last_k))
        }
    }
}

/// Returns the nominal segment propagation lengths `L_c` (`0..C`) of an `ems` `eme` problem on `mesh` without solving 2D modes.
pub fn eme_nominal_segment_lengths(
    mesh: &Mesh,
    spec: &nyx_physics::ProblemSpec,
    ctx: &nyx_physics::Ctx,
) -> Result<Vec<f64>, Error> {
    let problem = crate::problem::decode(spec, ctx)?;
    let Study::Eme { ref grid, y_range, n_modes, .. } = problem.study else {
        return Err(Error::invalid(format!(
            "eme_segment_length requires study 'eme', found '{}'",
            problem.study.kind()
        )));
    };
    let d3 = discretise(mesh)?;
    let (lo3, hi3) = d3.bbox();
    let [y0, y1] = y_range.unwrap_or([lo3[1], hi3[1]]);
    let mut segs = Vec::new();
    let mut sg = 0usize;
    let _ = build_grid_tree(grid, y0, y1, hi3[1], n_modes.max(1), &mut segs, &mut sg).map_err(Error::invalid)?;
    Ok(segs.into_iter().map(|s| s.length_m).collect())
}

/// Applies physical segment-length overrides `(segment_index, new_length_m)` onto `spec.study.params` of an `eme` `ProblemSpec`,
/// preserving the nominal 3D mesh sampling intervals (`mesh_lengths`) so cross-section cut planes never drift across adjacent regions.
pub fn eme_set_spec_segment_lengths(
    spec: &mut nyx_physics::ProblemSpec,
    ctx: &nyx_physics::Ctx,
    base_lengths: &[f64],
    overrides: &[(usize, f64)],
) -> Result<(), Error> {
    if overrides.is_empty() {
        return Ok(());
    }
    let problem = crate::problem::decode(spec, ctx)?;
    let Study::Eme { ref grid, .. } = problem.study else {
        return Err(Error::invalid(format!(
            "eme_segment_length requires study 'eme', found '{}'",
            problem.study.kind()
        )));
    };
    let mut ov_map = BTreeMap::new();
    for &(idx, val) in overrides {
        if !(val > 0.0 && val.is_finite()) {
            return Err(Error::invalid(format!("EME segment {idx} length must be positive and finite (got {val})")));
        }
        ov_map.insert(idx, val);
    }

    fn rewrite_grid(
        g: &EmeGridSpec,
        cursor: &mut usize,
        base_lengths: &[f64],
        ov: &BTreeMap<usize, f64>,
    ) -> (EmeGridSpec, bool) {
        match g {
            EmeGridSpec::Uniform { n_cells, n_modes, bend_radius } => {
                let nc = (*n_cells).max(1);
                let start = *cursor;
                *cursor += nc;
                let any = (start..start + nc).any(|i| ov.contains_key(&i));
                if !any {
                    return (g.clone(), false);
                }
                let m_lens: Vec<f64> =
                    (start..start + nc).map(|i| base_lengths.get(i).copied().unwrap_or(1e-6)).collect();
                let lens: Vec<f64> = (start..start + nc)
                    .map(|i| ov.get(&i).copied().unwrap_or_else(|| base_lengths.get(i).copied().unwrap_or(1e-6)))
                    .collect();
                (
                    EmeGridSpec::Explicit {
                        lengths: lens,
                        n_modes: n_modes.map(|m| vec![m; nc]),
                        bend_radius: bend_radius.map(|r| vec![Some(r); nc]),
                        mesh_lengths: Some(m_lens),
                    },
                    true,
                )
            }
            EmeGridSpec::Explicit { lengths, n_modes, bend_radius, mesh_lengths } => {
                let nc = lengths.len();
                let start = *cursor;
                *cursor += nc;
                let any = (start..start + nc).any(|i| ov.contains_key(&i));
                if !any {
                    return (g.clone(), false);
                }
                let m_lens = mesh_lengths.clone().unwrap_or_else(|| lengths.clone());
                let lens: Vec<f64> = (0..nc).map(|k| ov.get(&(start + k)).copied().unwrap_or(lengths[k])).collect();
                (
                    EmeGridSpec::Explicit {
                        lengths: lens,
                        n_modes: n_modes.clone(),
                        bend_radius: bend_radius.clone(),
                        mesh_lengths: Some(m_lens),
                    },
                    true,
                )
            }
            EmeGridSpec::Composite { subgrids } => {
                let mut changed = false;
                let mut out = Vec::with_capacity(subgrids.len());
                for sg in subgrids {
                    let (ng, ch) = rewrite_grid(sg, cursor, base_lengths, ov);
                    changed |= ch;
                    out.push(ng);
                }
                (EmeGridSpec::Composite { subgrids: out }, changed)
            }
            EmeGridSpec::PeriodicSubgrid { unit_cell, num_reps, unit_length } => {
                let (nuc, changed) = rewrite_grid(unit_cell, cursor, base_lengths, ov);
                let ul = if changed { None } else { *unit_length };
                (
                    EmeGridSpec::PeriodicSubgrid { unit_cell: Box::new(nuc), num_reps: *num_reps, unit_length: ul },
                    changed,
                )
            }
        }
    }

    let mut cur = 0usize;
    let (new_grid, _) = rewrite_grid(grid, &mut cur, base_lengths, &ov_map);
    validate_grid_spec(&new_grid).map_err(Error::invalid)?;
    let txt = serde_json::to_string(&new_grid).map_err(|e| Error::invalid(e.to_string()))?;
    if !spec.study.params.is_object() {
        spec.study.params = serde_json::json!({});
    }
    if let Some(obj) = spec.study.params.as_object_mut() {
        obj.insert("grid_kind".into(), serde_json::json!("json"));
        obj.insert("grid_json".into(), serde_json::json!(txt));
    }
    Ok(())
}

/// Evaluates [`nyx_physics::SensitivityReport`] for pure [`nyx_physics::DesignParam::EmeSegmentLength`] parameters
/// directly from cached [`EmeData`] in microseconds with zero 2D eigensolves.
pub fn evaluate_eme_length_sensitivity_from_data(
    data: &EmeData,
    objective: &nyx_physics::SensitivityObjective,
    params: &[nyx_physics::DesignParam],
    length_overrides: Option<&[f64]>,
) -> Result<nyx_physics::SensitivityReport, Error> {
    use nyx_physics::{DesignParam, SensitivityObjective, SensitivityReport};
    let zero_d = vec![0.0; data.total_segments];
    let (s_term, _) = data.evaluate_terminal_s_and_grad(length_overrides, &zero_d).map_err(Error::from_kernel)?;
    let p_ports = s_term.len();
    if p_ports == 0 {
        return Err(Error::invalid("eme sensitivity requires at least one modal or waveport channel"));
    }

    let mut aux_scalars = BTreeMap::new();
    let obj_val = match objective {
        SensitivityObjective::Transmission { out_port, in_port, target, .. } => {
            if *in_port >= p_ports || *out_port >= p_ports {
                return Err(Error::invalid(format!(
                    "transmission port index out of range (in_port={in_port}, out_port={out_port}, ports={p_ports})"
                )));
            }
            let mag_sq = s_term[*out_port][*in_port].norm_sqr();
            aux_scalars.insert("S_mag_sq".into(), mag_sq);
            match target {
                None => mag_sq,
                Some(t) => {
                    let diff = mag_sq - t;
                    0.5 * diff * diff
                }
            }
        }
        SensitivityObjective::SplitterBalance { in_port, out_ports, targets, .. } => {
            if *in_port >= p_ports || out_ports.iter().any(|&p| p >= p_ports) {
                return Err(Error::invalid("splitter_balance port index out of range"));
            }
            let m_out = out_ports.len();
            let eff_targets: Vec<f64> =
                if targets.is_empty() { vec![1.0 / m_out as f64; m_out] } else { targets.clone() };
            let mut obj = 0.0;
            for (mi, &p) in out_ports.iter().enumerate() {
                let mag_sq = s_term[p][*in_port].norm_sqr();
                aux_scalars.insert(format!("S_{p}_{in_port}_mag_sq"), mag_sq);
                obj += 0.5 * (mag_sq - eff_targets[mi]).powi(2);
            }
            obj
        }
        SensitivityObjective::EffectiveIndex { mode, target_neff } => {
            let us0 = data.unique_slices.first().ok_or_else(|| Error::invalid("EME cascade has zero slices"))?;
            let &[re, im] = us0.neff_fwd.get(*mode).ok_or_else(|| {
                Error::invalid(format!("requested mode {mode}, EME slice solved {} mode(s)", us0.neff_fwd.len()))
            })?;
            aux_scalars.insert("neff_re".into(), re);
            aux_scalars.insert("neff_im".into(), im);
            match target_neff {
                None => re,
                Some(nt) => 0.5 * (re - nt).powi(2),
            }
        }
        other => {
            return Err(Error::Unsupported(format!(
                "ems eme sensitivity supports Transmission, SplitterBalance, and EffectiveIndex, got '{other:?}'"
            )));
        }
    };

    let contract_ds = |ds_term: &[Vec<C>]| -> f64 {
        match objective {
            SensitivityObjective::Transmission { out_port, in_port, target, .. } => {
                let s_pq = s_term[*out_port][*in_port];
                let d_mag_sq = 2.0 * (s_pq.conj() * ds_term[*out_port][*in_port]).re;
                match target {
                    None => d_mag_sq,
                    Some(t) => (s_pq.norm_sqr() - t) * d_mag_sq,
                }
            }
            SensitivityObjective::SplitterBalance { in_port, out_ports, targets, .. } => {
                let m_out = out_ports.len();
                let eff_targets: Vec<f64> =
                    if targets.is_empty() { vec![1.0 / m_out as f64; m_out] } else { targets.clone() };
                let mut g = 0.0;
                for (mi, &p) in out_ports.iter().enumerate() {
                    let s_pq = s_term[p][*in_port];
                    let d_mag_sq = 2.0 * (s_pq.conj() * ds_term[p][*in_port]).re;
                    g += (s_pq.norm_sqr() - eff_targets[mi]) * d_mag_sq;
                }
                g
            }
            _ => 0.0,
        }
    };

    let mut gradients = Vec::with_capacity(params.len());
    let mut flat_gradient = Vec::with_capacity(params.len());
    let mut param_names = Vec::with_capacity(params.len());

    for param in params {
        param_names.push(param.label());
        let DesignParam::EmeSegmentLength { segment, .. } = param else {
            return Err(Error::invalid(
                "evaluate_eme_length_sensitivity_from_data expects only EmeSegmentLength params",
            ));
        };
        if *segment >= data.total_segments {
            return Err(Error::invalid(format!(
                "eme_segment_length segment {segment} out of range (cascade has {} segment(s))",
                data.total_segments
            )));
        }
        let g = if matches!(objective, SensitivityObjective::EffectiveIndex { .. }) {
            0.0
        } else {
            let mut d_lengths = vec![0.0; data.total_segments];
            d_lengths[*segment] = 1.0;
            let (_, ds_term) =
                data.evaluate_terminal_s_and_grad(length_overrides, &d_lengths).map_err(Error::from_kernel)?;
            contract_ds(&ds_term)
        };
        gradients.push(vec![g]);
        flat_gradient.push(g);
    }

    let rep = SensitivityReport {
        objective: obj_val,
        gradients,
        flat_gradient,
        param_names,
        aux_scalars,
        projected_densities: BTreeMap::new(),
        volume_fractions: BTreeMap::new(),
        volume_fraction_gradients: BTreeMap::new(),
        mode_basis: None,
        log: Vec::new(),
    };
    rep.check_finite()?;
    Ok(rep)
}

// ------------------------------------------------------------------------------------------------------------- Stage 3: Fields & SolveResult

/// Stage 3: packages [`DrivenPortData`], [`EmeData`], solver diagnostics, and the reconstructed 3D electric field
/// into a [`SolveResult`] (`ADR 9`, `ADR 140`).
pub fn eme_fields(_mesh: &Mesh, asm: &EmeAssembly, cascade: &EmeCascade) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let d3 = &asm.mesh_3d;
    let order = asm.problem.order;
    let p_count = cascade.data.ports.len();

    // Build DrivenPortData (single frequency `[freq_hz]`, compatible with `nyx_analysis::network` and Touchstone `.sNp`)
    let (s_pairs, z_pairs, y_pairs, conv_rep) = s_to_z_and_y_checked(&cascade.s_terminal, &cascade.data.port_z0_ohm);
    let (passivity_max_sv, reciprocity_err) = s_matrix_passivity_and_reciprocity(&cascade.s_terminal);

    let driven_data = DrivenPortData {
        freqs_hz: vec![asm.freq_hz],
        ports: cascade.data.ports.clone(),
        z0_ohm: vec![cascade.data.port_z0_ohm.iter().map(|&z| [z, 0.0]).collect()],
        s_matrix: vec![s_pairs],
        z_matrix: vec![z_pairs],
        y_matrix: vec![y_pairs],
        hdm_solves: asm.unique_slices.len(),
        rom_dim: asm.unique_slices.iter().map(|u| u.max_modes).sum(),
        anchor_freqs_hz: vec![asm.freq_hz],
        max_residual: cascade.max_residual,
    };

    let (fields, cell_fields) = if asm.reconstruct_fields {
        reconstruct_3d_eme_fields(asm, cascade)?
    } else {
        let nd = d3.n_dofs(order);
        let z = || vec![0.0; nd];
        (
            vec![
                Field::scalar("|E|", "a.u.", z()),
                Field::scalar("Re Ex", "a.u.", z()),
                Field::scalar("Re Ey", "a.u.", z()),
                Field::scalar("Re Ez", "a.u.", z()),
            ],
            Vec::new(),
        )
    };

    let mut scalars = vec![
        ("unique_slices".into(), Scalar::Num(asm.unique_slices.len() as f64)),
        ("total_segments".into(), Scalar::Num(asm.mesh_segments.len() as f64)),
        ("passivity_max_sigma".into(), Scalar::Num(passivity_max_sv)),
        ("reciprocity_rel_err".into(), Scalar::Num(reciprocity_err)),
    ];
    if p_count >= 1 {
        scalars.push(("R_p1".into(), Scalar::Num(cascade.s_terminal[0][0].norm_sqr())));
    }
    if p_count >= 2 {
        let s21 = cascade.s_terminal[1][0];
        scalars.push(("T_21".into(), Scalar::Num(s21.norm_sqr())));
        scalars.push(("phase_S21_rad".into(), Scalar::Num(s21.im.atan2(s21.re))));
    }
    for (u_idx, sm) in cascade.slice_modes.iter().enumerate() {
        if let Some(nf0) = sm.neff_fwd.first() {
            scalars.push((format!("slice_{}_neff1_re", u_idx + 1), Scalar::Num(nf0.re)));
            scalars.push((format!("slice_{}_neff1_im", u_idx + 1), Scalar::Num(nf0.im)));
        }
    }
    scalars.push(crate::optical::mode_frame(&asm.problem).scalar());
    scalars.push(driven_data.scalar());
    scalars.push(cascade.data.scalar()?);

    let mut log = vec![format!(
        "EME cascade at f = {:.6e} Hz: {} segment(s), {} unique slice(s) (BLAKE3 deduplicated), {} terminal port(s), sigma_max(S) = {:.6}",
        asm.freq_hz,
        asm.mesh_segments.len(),
        asm.unique_slices.len(),
        p_count,
        passivity_max_sv
    )];
    for w in conv_rep.warning_messages {
        log.push(format!("network warning: {w}"));
    }

    let mut diag = crate::diag::EmsDiag::new("eme_redheffer_cascade");
    let all_res: Vec<f64> = cascade.slice_modes.iter().flat_map(|s| s.residuals.iter().copied()).collect();
    diag.record_mode_residuals(
        "eme_redheffer_cascade",
        if all_res.is_empty() { vec![0.0] } else { all_res },
        cascade.all_converged,
        Some(cascade.total_restarts),
    );
    let (diagnostics, error_budget) = diag.finish();

    let sol = Solution {
        label: format!("EME f = {:.4e} Hz (port 1 drive)", asm.freq_hz),
        fields,
        cell_fields,
        value: Some(asm.freq_hz),
        eigenvalue: None,
    };

    Ok(SolveResult {
        physics: "ems".into(),
        study: "eme".into(),
        order: order as u8,
        dim: 3,
        solutions: vec![sol],
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: asm.slice_modes_total_dofs(),
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}

impl EmeAssembly {
    fn slice_modes_total_dofs(&self) -> usize {
        self.unique_slices.iter().map(|u| u.mesh_2d.n_dofs(self.problem.order)).sum()
    }
}

/// Reconstructs the bounded 3D electric field `E(x, y, z)` across `asm.mesh_3d` under unit incident drive on port 1.
fn reconstruct_3d_eme_fields(
    asm: &EmeAssembly,
    cascade: &EmeCascade,
) -> Result<(Vec<Field>, Vec<nyx_physics::CellField>), String> {
    let d3 = &asm.mesh_3d;
    let order = asm.problem.order;
    let nc3 = d3.n_cells();
    let nloc3 = d3.nloc(order);
    let n_segs = asm.mesh_segments.len();

    // Compute forward a_c(y_lo) and backward b_c(y_hi) modal amplitudes on each mesh segment via prefix/suffix Redheffer blocks
    let (a_seg, b_seg) = solve_mesh_segment_amplitudes(asm, cascade)?;

    // Build TransferIndex + SliceCellEval per unique 2D slice
    let empty_names = BTreeMap::new();
    let opts = TransferOptions { snap_rel: 1e-3, parallel: false, ..TransferOptions::default() };
    let mut slice_indices = Vec::with_capacity(cascade.slice_modes.len());
    let mut slice_evals = Vec::with_capacity(cascade.slice_modes.len());
    for sm in &cascade.slice_modes {
        slice_indices
            .push(TransferIndex::new(SourceMesh::Simplex(&sm.mesh_2d), &empty_names, opts).map_err(|e| e.to_string())?);
        slice_evals.push(SliceCellEval::new(sm));
    }

    let bary3 = local_bary(3, order);
    let mut ex_c = vec![c(0.0, 0.0); nc3 * nloc3];
    let mut ey_c = vec![c(0.0, 0.0); nc3 * nloc3];
    let mut ez_c = vec![c(0.0, 0.0); nc3 * nloc3];
    let mut mode_scratch = Vec::new();
    let mut pd2 = Vec::new();

    for ci in 0..nc3 {
        let cv = d3.cell(ci);
        let pts =
            [d3.point(cv[0] as usize), d3.point(cv[1] as usize), d3.point(cv[2] as usize), d3.point(cv[3] as usize)];
        let yc = 0.25 * (pts[0][1] + pts[1][1] + pts[2][1] + pts[3][1]);
        // Find owning mesh segment by y-centroid
        let seg_i = find_segment_for_y(&asm.mesh_segments, yc);
        let seg = &asm.mesh_segments[seg_i];
        let u_idx = seg.unique_slice;
        let sm = &cascade.slice_modes[u_idx];
        let idx2 = &slice_indices[u_idx];
        let eval2 = &slice_evals[u_idx];
        let k_m = a_seg[seg_i].len().min(sm.neff_fwd.len());
        let span_seg = (seg.y_mesh_hi - seg.y_mesh_lo).max(1e-30);
        let y_scale = seg.length_m / span_seg;

        for (a_loc, lam3) in bary3.iter().enumerate() {
            let x = lam3[0] * pts[0][0] + lam3[1] * pts[1][0] + lam3[2] * pts[2][0] + lam3[3] * pts[3][0];
            let y = lam3[0] * pts[0][1] + lam3[1] * pts[1][1] + lam3[2] * pts[2][1] + lam3[3] * pts[3][1];
            let z = lam3[0] * pts[0][2] + lam3[1] * pts[1][2] + lam3[2] * pts[2][2] + lam3[3] * pts[3][2];

            let Some(loc2) = idx2.locate(&[x, z], None) else { continue };
            let lam2 = [loc2.coords[0], loc2.coords[1], loc2.coords[2]];
            eval2.eval_all_modes(loc2.cell, &lam2, &mut mode_scratch);
            let (nv2, _) = shape(&lam2, order);
            sm.mesh_2d.cell_dofs(loc2.cell, order, &mut pd2);

            let dy_fwd = ((y - seg.y_mesh_lo) * y_scale).clamp(0.0, seg.length_m);
            let dy_bwd = ((seg.y_mesh_hi - y) * y_scale).clamp(0.0, seg.length_m);

            let (mut ex, mut ey, mut ez) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            for m in 0..k_m {
                let pf = a_seg[seg_i][m] * stable_prop_factor([sm.beta_fwd[m].re, sm.beta_fwd[m].im], dy_fwd);
                let pb = b_seg[seg_i][m] * stable_prop_factor([sm.beta_bwd[m].re, sm.beta_bwd[m].im], dy_bwd);
                let amp_t = pf + pb;
                let amp_y = pf - pb;
                ex += amp_t * mode_scratch[m][0];
                ez += amp_t * mode_scratch[m][1];

                let mut phi_q = c(0.0, 0.0);
                for (k_n, &dof) in pd2.iter().enumerate() {
                    phi_q += sm.phi[m][dof] * nv2[k_n];
                }
                let ey_m = c(0.0, 1.0) * sm.beta_fwd[m] * phi_q;
                ey += amp_y * ey_m;
            }
            let flat = ci * nloc3 + a_loc;
            ex_c[flat] = ex;
            ey_c[flat] = ey;
            ez_c[flat] = ez;
        }
    }

    let _ = n_segs;
    let big = ex_c.iter().chain(ey_c.iter()).chain(ez_c.iter()).copied().fold(c(0.0, 0.0), |m, z| {
        if z.norm() > m.norm() {
            z
        } else {
            m
        }
    });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };

    let mag_cell: Vec<f64> =
        (0..nc3 * nloc3).map(|i| (ex_c[i].norm_sqr() + ey_c[i].norm_sqr() + ez_c[i].norm_sqr()).sqrt()).collect();
    let mag_nodal = node_average(d3, &mag_cell, 1, order);
    let scale = 1.0 / mag_nodal.iter().copied().fold(0.0_f64, f64::max).max(1e-300);

    let re_scaled = |arr: &[C]| -> Vec<f64> {
        let r_cell: Vec<f64> = arr.iter().map(|&z| (z * rot).re * scale).collect();
        node_average(d3, &r_cell, 1, order)
    };
    let mut cells = CellFields::new(d3, order);
    for (name, arr) in [("Ex", &ex_c), ("Ey", &ey_c), ("Ez", &ez_c)] {
        cells.push_complex(
            name,
            "a.u.",
            arr.iter().map(|&z| (z * rot).re * scale).collect(),
            arr.iter().map(|&z| (z * rot).im * scale).collect(),
        );
    }

    let fields = vec![
        Field::scalar("|E|", "a.u.", mag_nodal.into_iter().map(|v| v * scale).collect()),
        Field::scalar("Re Ex", "a.u.", re_scaled(&ex_c)),
        Field::scalar("Re Ey", "a.u.", re_scaled(&ey_c)),
        Field::scalar("Re Ez", "a.u.", re_scaled(&ez_c)),
    ];
    Ok((fields, cell_fields_of(&cells)))
}

fn find_segment_for_y(segs: &[EmeSegmentPlan], y: f64) -> usize {
    for (i, s) in segs.iter().enumerate() {
        if y <= s.y_mesh_hi || i + 1 == segs.len() {
            return i;
        }
    }
    0
}

/// Computes bounded internal modal wave amplitudes `(a_c(0), b_c(L_c))` on each mesh segment `0..C`
/// under unit excitation of the first input port channel.
fn solve_mesh_segment_amplitudes(
    asm: &EmeAssembly,
    cascade: &EmeCascade,
) -> Result<(Vec<Vec<C>>, Vec<Vec<C>>), String> {
    let segs = &asm.mesh_segments;
    let c_len = segs.len();
    let data = &cascade.data;

    let mut if_lookup: HashMap<(usize, usize), ScatterBlock> = HashMap::new();
    for ov in &data.interfaces {
        if_lookup.insert(
            (ov.left_id, ov.right_id),
            ScatterBlock { n_left: ov.n_left, n_right: ov.n_right, s: pairs_to_c_mat(&ov.s_interface) },
        );
    }

    let seg_blocks: Vec<ScatterBlock> = segs
        .iter()
        .map(|s| {
            let us = &data.unique_slices[s.unique_slice];
            let k = s.n_modes.min(us.n_modes).max(1);
            propagation_block(&us.beta_fwd, &us.beta_bwd, k, s.length_m)
        })
        .collect();

    // Ensure adjacent mesh-segment interfaces exist (e.g. if mesh_segments came from a single Periodic unit cell, all adjacent pairs are already in `data.interfaces`)
    let mut adj_if = Vec::with_capacity(c_len.saturating_sub(1));
    for i in 0..c_len.saturating_sub(1) {
        let (ul, ur) = (segs[i].unique_slice, segs[i + 1].unique_slice);
        let (kl, kr) = (seg_blocks[i].n_right, seg_blocks[i + 1].n_left);
        if let Some(blk) = if_lookup.get(&(ul, ur)) {
            if blk.n_left == kl && blk.n_right == kr {
                adj_if.push(blk.clone());
                continue;
            }
        }
        let ov = compute_interface_overlap(ul, &cascade.slice_modes[ul], ur, &cascade.slice_modes[ur])?;
        let s_if = build_interface_s_matrix(&ov.d_left, &ov.d_right_tilde, &ov.v_lr, &ov.w_lr, kl, kr)?;
        adj_if.push(ScatterBlock { n_left: kl, n_right: kr, s: s_if });
    }

    // Prefix blocks `pre_left[c]` = cascade from external input up to left end of segment `c`
    let m0 = seg_blocks[0].n_left;
    let init_pre = if let Some(ref in_ov) = data.input_port_overlap {
        ScatterBlock { n_left: in_ov.n_left, n_right: in_ov.n_right, s: pairs_to_c_mat(&in_ov.s_interface) }
    } else {
        identity_scatter_block(m0)
    };

    let mut pre_left = Vec::with_capacity(c_len);
    let mut cur_pre = init_pre;
    for i in 0..c_len {
        pre_left.push(cur_pre.clone());
        if i + 1 < c_len {
            let after_seg = redheffer_star(&cur_pre, &seg_blocks[i])?;
            cur_pre = redheffer_star(&after_seg, &adj_if[i])?;
        }
    }

    // Suffix blocks `suf_right[c]` = cascade from right end of segment `c` to external output
    let m_last = seg_blocks[c_len - 1].n_right;
    let init_suf = if let Some(ref out_ov) = data.output_port_overlap {
        ScatterBlock { n_left: out_ov.n_left, n_right: out_ov.n_right, s: pairs_to_c_mat(&out_ov.s_interface) }
    } else {
        identity_scatter_block(m_last)
    };

    let mut suf_right = vec![init_suf.clone(); c_len];
    let mut cur_suf = init_suf;
    for i in (0..c_len.saturating_sub(1)).rev() {
        let with_next = redheffer_star(&seg_blocks[i + 1], &cur_suf)?;
        cur_suf = redheffer_star(&adj_if[i], &with_next)?;
        suf_right[i] = cur_suf.clone();
    }

    // Determine incident input vector `a_ext` (port 0 = 1, others = 0)
    let n_ext_in = pre_left[0].n_left;
    let mut a_ext = vec![c(0.0, 0.0); n_ext_in];
    if let Some(&(is_out, ch)) = data.port_channels.first() {
        if !is_out && ch < n_ext_in {
            a_ext[ch] = c(1.0, 0.0);
        } else if n_ext_in > 0 {
            a_ext[0] = c(1.0, 0.0);
        }
    } else if n_ext_in > 0 {
        a_ext[0] = c(1.0, 0.0);
    }

    let mut a_seg = Vec::with_capacity(c_len);
    let mut b_seg = Vec::with_capacity(c_len);
    for i in 0..c_len {
        let k = seg_blocks[i].n_left;
        let pre = &pre_left[i];
        let suf_from_left = redheffer_star(&seg_blocks[i], &suf_right[i])?;
        // (I - S22^pre * S11^suf_left) a_i(0) = S21^pre * a_ext
        let mut mat = vec![c(0.0, 0.0); k * k];
        for r in 0..k {
            mat[r * k + r] = c(1.0, 0.0);
            for col in 0..k {
                let mut sum = c(0.0, 0.0);
                for q in 0..k {
                    sum += pre.s[pre.n_left + r][pre.n_left + q] * suf_from_left.s[q][col];
                }
                mat[r * k + col] -= sum;
            }
        }
        let mut rhs = vec![c(0.0, 0.0); k];
        for r in 0..k {
            for col in 0..pre.n_left {
                rhs[r] += pre.s[pre.n_left + r][col] * a_ext[col];
            }
        }
        let a_i = solve_with_reg(&mat, &[rhs], k)?.remove(0);

        // At right end of segment i: forward arriving wave is a_i(L_i) = P_fwd * a_i(0), and backward entering wave is b_i(L_i) = S11^{suf_right[i]} * a_i(L_i)
        let us = &data.unique_slices[segs[i].unique_slice];
        let mut a_right = vec![c(0.0, 0.0); k];
        for m in 0..k {
            a_right[m] = a_i[m] * stable_prop_factor(us.beta_fwd[m], segs[i].length_m);
        }
        let mut b_i_right = vec![c(0.0, 0.0); k];
        for r in 0..k {
            for col in 0..k {
                b_i_right[r] += suf_right[i].s[r][col] * a_right[col];
            }
        }
        a_seg.push(a_i);
        b_seg.push(b_i_right);
    }

    Ok((a_seg, b_seg))
}

fn identity_scatter_block(m: usize) -> ScatterBlock {
    let n = 2 * m;
    let mut s = vec![vec![c(0.0, 0.0); n]; n];
    for i in 0..m {
        s[m + i][i] = c(1.0, 0.0);
        s[i][m + i] = c(1.0, 0.0);
    }
    ScatterBlock { n_left: m, n_right: m, s }
}

/// One-call entry for `Study::Eme`: runs Stages 1, 2, and 3 in sequence (`ADR 9`, `ADR 140`).
pub fn solve_eme(mesh: &Mesh, problem: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let asm = assemble_eme(mesh, problem)?;
    femcore::progress::check()?;
    let cascade = solve_eme_cascade(&asm)?;
    femcore::progress::check()?;
    let mut res = eme_fields(mesh, &asm, &cascade)?;
    res.seconds = t0.elapsed().as_secs_f64();
    Ok(res)
}
